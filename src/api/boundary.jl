###############################################################
# The plain-array API boundary
###############################################################
#
# Pagos' public API deals in plain arrays and scalars; Chmy `Field`s are an internal
# representation (see `roadmaps/chmy.md`, Phase 2). Two functions carry data across that
# boundary in either direction, and both work on *either* representation, so calling code
# does not need to know which one a state struct was built with:
#
#   asarray(x)          Field or array  ->  plain array (a view, no copy)
#   setdata!(dst, src)  array or scalar ->  Field or array (a checked copy)
#
# Neither takes a grid, a `Runtime` or any Chmy type, and neither is specific to a node
# class: they move data, they do not interpret staggering.

"""
$(TYPEDSIGNATURES)

The data of `x` as a plain `AbstractArray`, whichever representation `x` is in: the
interior view for a Chmy `Field` (halos excluded), and `x` itself for an array. This is
the accessor the public API and any output path (I/O, plotting, diagnostics) should use —
it is a `view`, not a copy, so writes through it land in the field.

Use it in preference to broadcasting a `Field` directly. A Chmy `Field` *is* an
`AbstractArray` (`Chmy.AbstractField{T,N,L} <: AbstractArray{T,N}`), so `f .+ 1` compiles,
but it has no `BroadcastStyle` of its own and so falls back to generic scalar `getindex` —
correct on CPU, catastrophically slow or broken on a GPU. `asarray(f)` yields a view of
the real underlying array, which broadcasts natively on every backend.

!!! note "Lazy fields are returned as-is"
    For a `Chmy.FunctionField` or `ConstantField` there is no stored array, and Chmy's own
    `interior` is the identity — so `asarray` returns the lazy field unchanged. Every field
    in a Pagos state struct is an allocated `Chmy.Field`, so this only arises for fields a
    caller constructed themselves; `collect` them first if a materialized array is needed.

# Examples

```jldoctest
julia> grid = StaggeredGrid(Float64, 4.0, 4.0, 1.0, 1.0);

julia> topo = TopographicState(grid);

julia> size(asarray(topo.thickness.ice))
(4, 4, 1)

julia> asarray(zeros(3, 3)) == zeros(3, 3)      # a no-op on plain arrays
true
```
"""
asarray(f::Chmy.AbstractField) = interior(f)
asarray(A::AbstractArray) = A

# Trailing singleton dimensions carry no data, so a depth-integrated field (`nx, ny, 1`)
# and the `nx, ny` matrix a user naturally has for it describe the same thing. Comparing
# sizes with those dimensions dropped accepts that pairing while still rejecting a genuine
# mismatch.
function _data_size(sz::Tuple)
    n = length(sz)
    while n > 1 && sz[n] == 1
        n -= 1
    end
    return ntuple(i -> sz[i], n)
end

function _check_data_shape(dst_size, src_size)
    _data_size(dst_size) == _data_size(src_size) && return nothing
    throw(
        DimensionMismatch(
            "cannot write data of size $src_size into a destination of size $dst_size. " *
            "Trailing singleton dimensions are ignored, so e.g. an (nx, ny) matrix is " *
            "accepted for an (nx, ny, 1) depth-integrated field, but the leading " *
            "dimensions must agree exactly. If the source was dimensioned against a " *
            "`RegularGrid`, note that `StaggeredGrid` is one cell smaller per horizontal " *
            "axis by construction (see the `StaggeredGrid` docstring) — re-dimension the " *
            "data rather than padding it.",
        ),
    )
end

"""
$(TYPEDSIGNATURES)

Write `src` — a plain array, a scalar, or another field — into `dst`, which may be a Chmy
`Field` or a plain array. Only `dst`'s interior is touched; halos are left alone (they are
filled by boundary conditions, see `roadmaps/chmy.md`, Phase 4). Element types are
converted as needed, so a `Float64` array initializes a `Float32` field, and a host array
initializes a device field by copy. Returns `dst`.

This is the promotion half of the plain-array boundary, and it exists because
**`Chmy.set!(f, A::AbstractArray)` is an unchecked `copyto!`**: it copies by *linear*
index and only errors when `A` is too large. An array that is too small silently fills
part of the field and leaves the rest at whatever was there before; one of the right
length but the wrong shape is silently reinterpreted. Both are easy to hit at exactly this
boundary — `StaggeredGrid` is deliberately one cell smaller per horizontal axis than
`RegularGrid` for the same `(lx, ly, dx, dy)`, so data carried over from the old grid is
off by one in each direction. `setdata!` checks the shape first and raises
`DimensionMismatch` naming both sizes.

# Examples

```jldoctest
julia> grid = StaggeredGrid(Float64, 4.0, 4.0, 1.0, 1.0);

julia> topo = TopographicState(grid);

julia> setdata!(topo.thickness.ice, 1000.0);        # a scalar

julia> setdata!(topo.thickness.ice, fill(2000.0, 4, 4));   # an (nx, ny) matrix

julia> asarray(topo.thickness.ice)[1, 1, 1]
2000.0

julia> topo.thickness.ice isa Chmy.Fields.Field    # the caller never had to name this
true
```

A source of any other shape — `fill(1.0, 5, 5)` here, or anything sized against a
`RegularGrid` — raises `DimensionMismatch` naming both sizes rather than being copied.
"""
function setdata!(dst, src::AbstractArray)
    d = asarray(dst)
    _check_data_shape(size(d), size(src))
    copyto!(d, size(d) == size(src) ? src : reshape(src, size(d)))
    return dst
end

function setdata!(dst, src::Number)
    fill!(asarray(dst), src)
    return dst
end

# Routed through `asarray` rather than the method above so the read side is a plain view
# too: `copyto!` from a `Field` would go through its scalar `getindex`, which is the same
# GPU trap `asarray` exists to avoid.
setdata!(dst, src::Chmy.AbstractField) = setdata!(dst, asarray(src))

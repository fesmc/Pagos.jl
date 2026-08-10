###############################################################
# Device-native bulk fill from a host source array
###############################################################
#
# The naive way to seed a `Field` from external data (a NetCDF restart, say) is an
# element-by-element loop: `f[i, j, k] = data[i, j]` over the interior, plus a clamped
# read for the halo ring so the field still has sensible values just outside the domain.
# That loop needs scalar `getindex`/`setindex!` on `f`, which is not available on a
# `CuArray`-backed field — so the field has to be built on the host, filled there, and
# only then transferred to the device, doubling peak memory for the whole state (see
# `pagos-roadmap/memreduce.md`, "Build device states without a host duplicate").
#
# `copyto!(interior(f), data)` has none of that problem: it works on a device field
# without scalar indexing (it lowers to a single H→D memcpy). The only thing missing is
# an equally bulk way to fill the halo ring with the same "clamp to the nearest interior
# column" extension the scalar loop used — which is what this file adds, by building the
# clamped extension as an ordinary (small, host-side) array operation and handing the
# whole padded block to `copyto!` in one call.

_as_ntuple(h::Tuple, ::Integer) = h
_as_ntuple(h::Integer, n::Integer) = ntuple(_ -> h, n)

"""
$(TYPEDSIGNATURES)

The number of ghost cells actually reachable on `f`'s first two (x, y) axes, as `(hx,
hy)` — regardless of whether `halo(f)` is a single `Int` (uniform on every axis) or a
per-axis tuple (as produced by [`_halo2d`](@ref) for depth-integrated fields).

Chmy's `halo` type parameter is *not* itself a ring count: `Field` allocates
`size(grid, loc) .+ 4 .* halo` (see `Chmy.Fields.Field`), i.e. `2 * halo` ghost rings per
side. `2 .* halo(f)` is therefore the actual x/y padding width to fill.
"""
_xy_halo(f::Chmy.AbstractField) = 2 .* _as_ntuple(halo(f), 2)[1:2]

"""
$(TYPEDSIGNATURES)

A view into `f`'s backing array covering the interior expanded by `hx`/`hy` on the first
two axes and left at the bare interior on every other axis (in particular: the z axis, if
any, is never touched). Mirrors Chmy's own `interior(f; with_halo=true)`, but expands only
the axes the caller asks for, so a field's z halo (if it has one) is left exactly as
allocated rather than being read as if it were meaningful padding.
"""
function _xy_padded_interior(f::Chmy.AbstractField, hx::Integer, hy::Integer)
    ax = axes(f)
    Ht = _as_ntuple(halo(f), length(ax))
    rngs = ntuple(length(ax)) do d
        if d == 1
            (first(ax[1])-hx):(last(ax[1])+hx)
        elseif d == 2
            (first(ax[2])-hy):(last(ax[2])+hy)
        else
            ax[d]
        end
    end
    idx = ntuple(d -> rngs[d] .+ 2 * Ht[d], length(ax))
    return view(parent(f), idx...)
end

# Clamped x/y extension of a small host array, e.g. `data[ic, jc]` with `ic = clamp(i, 1,
# ni)` — computed as ordinary fancy indexing rather than a loop, since `data` is host-side
# and tiny (tens of MB) regardless of where `f` itself lives.
function _clamp_pad_xy(data::AbstractArray, hx::Integer, hy::Integer)
    ni, nj = size(data, 1), size(data, 2)
    ix = clamp.((1-hx):(ni+hx), 1, ni)
    iy = clamp.((1-hy):(nj+hy), 1, nj)
    return data[ix, iy, ntuple(_ -> Colon(), ndims(data) - 2)...]
end

"""
$(TYPEDSIGNATURES)

Fill `f`'s interior from the 2D array `data` (broadcasting it over every z level) and its
x/y halo ring by clamping to the nearest interior column/row — the same extension a naive
`f[i, j, k] = data[clamp(i, 1, ni), clamp(j, 1, nj)]` loop over the halo would produce, but
built as two bulk array operations (`data[ix, iy]` then `copyto!`) so it works on a
`CuArray`-backed `f` without ever touching it with scalar indexing. `f`'s z halo, if any,
is left untouched, matching the scalar version (which never wrote it either).

Use this in place of a hand-rolled element-by-element fill when seeding state directly on
the device — see `pagos-roadmap/memreduce.md`, "Build device states without a host
duplicate".
"""
function fill_from_grid!(f::Chmy.AbstractField, data::AbstractMatrix)
    hx, hy = _xy_halo(f)
    padded = _clamp_pad_xy(data, hx, hy)
    dst = _xy_padded_interior(f, hx, hy)
    copyto!(dst, repeat(padded, 1, 1, size(dst, 3)))
    return f
end

"""
$(TYPEDSIGNATURES)

Three-dimensional counterpart of [`fill_from_grid!`](@ref): `data` already carries its own
z levels (shape `(ni, nj, nk)`, matching `f`'s interior exactly, no broadcasting over z),
so only the x/y halo ring needs the clamped extension.
"""
function fill_from_grid3d!(f::Chmy.AbstractField, data::AbstractArray{<:Any,3})
    hx, hy = _xy_halo(f)
    padded = _clamp_pad_xy(data, hx, hy)
    dst = _xy_padded_interior(f, hx, hy)
    copyto!(dst, padded)
    return f
end

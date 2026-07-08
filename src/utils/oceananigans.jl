"""
    StaggeredGrids

A tiny, self-contained re-implementation of the slice of Oceananigans we actually
need: a staggered `RectilinearGrid`, a located `Field`, the finite-difference /
interpolation operators, and the `@at` macro.

The whole thing rests on two ideas:

  1. **Every node carries a location** `(LX, LY, LZ)` with each entry `Center`,
     `Face` (or `Nothing` for a reduced axis). Leaves (`Field`s) get it from
     construction; operators compute theirs from their operands.
  2. **Every node is index-evaluable**: `evaluate(node, i, j, k, grid)` returns the
     value at grid point `(i, j, k)` *at that node's own location*. A `Field` reads
     its array (with topology-aware index wrapping); an operator recurses into its
     operands at shifted indices. Interpolation is just another node, inserted
     automatically whenever an operator needs an operand at a location it isn't at.

# Example

```julia
using .StaggeredGrids

grid = RectilinearGrid(topology = (Periodic, Flat, Bounded),
                       size = (4, 4),
                       x = (0, 2π),
                       z = (-4, 0))

v_x = Field{Face, Center, Center}(grid)
v_x.data .= rand(size(v_x.data)...)

op  = @at (Center, Center, Center) ∂x(v_x)^2 + ∂y(v_x)^2
out = compute(op)          # a Field at (Center, Center, Center)
```

This is a first, CPU-only cut: boundaries are handled by topology-aware index
wrapping inside `evaluate(::Field, ...)` (periodic wrap / bounded clamp / flat →1),
so no halo bookkeeping is needed yet. Spacing is uniform. See the TODOs at the
bottom for the natural extensions (stretched grids, halos, KernelAbstractions).
"""
module StaggeredGrids

import Adapt

export RectilinearGrid, Field
export Center, Face, Periodic, Bounded, Flat
export ∂x, ∂y, ∂z, @at, compute, compute!, location, evaluate, interpolate

# ---------------------------------------------------------------------------
# 1. Locations and topologies (singleton *types*, used as type parameters)
# ---------------------------------------------------------------------------

abstract type AbstractLocation end
struct Center <: AbstractLocation end
struct Face   <: AbstractLocation end
# `Nothing` doubles as the location/topology of a reduced (Flat) axis.

abstract type AbstractTopology end
struct Periodic <: AbstractTopology end
struct Bounded  <: AbstractTopology end
struct Flat     <: AbstractTopology end

flip(::Type{Center})  = Face
flip(::Type{Face})    = Center
flip(::Type{Nothing}) = Nothing

# Number of grid points carried by a (topology, location) pair along an axis of
# `N` cells. This is the *only* place the staggering shows up in array sizes.
npoints(::Type{Periodic}, ::Type{Center}, N) = N
npoints(::Type{Periodic}, ::Type{Face},   N) = N        # Nth+1 face ≡ 1st
npoints(::Type{Bounded},  ::Type{Center}, N) = N
npoints(::Type{Bounded},  ::Type{Face},   N) = N + 1    # extra boundary face
npoints(::Type{Flat},     _,              N) = 1

# ---------------------------------------------------------------------------
# 2. RectilinearGrid (uniform spacing for now)
# ---------------------------------------------------------------------------

struct RectilinearGrid{TX, TY, TZ, FT}
    Nx::Int; Ny::Int; Nz::Int
    Δx::FT;  Δy::FT;  Δz::FT
end

topologies(::RectilinearGrid{TX, TY, TZ}) where {TX, TY, TZ} = (TX, TY, TZ)
Base.eltype(::RectilinearGrid{TX, TY, TZ, FT}) where {TX, TY, TZ, FT} = FT

"""
    RectilinearGrid([FT=Float64]; topology, size, x=nothing, y=nothing, z=nothing)

`topology` is a 3-tuple of `Periodic`/`Bounded`/`Flat`. `size` lists the number of
cells of the **non-Flat** dimensions, in order (Oceananigans convention). `x`/`y`/`z`
give the extent `(lo, hi)` of each non-Flat dimension; spacing is `extent / N`.
"""
function RectilinearGrid(FT::Type{<:AbstractFloat} = Float64;
                         topology, size, x = nothing, y = nothing, z = nothing)
    TX, TY, TZ = topology
    tops  = (TX, TY, TZ)
    sizes = Tuple(size)

    # Distribute the `size` tuple across the non-Flat dimensions, in order.
    N = ones(Int, 3)
    c = 1
    for d in 1:3
        if tops[d] !== Flat
            N[d] = sizes[c]
            c += 1
        end
    end
    Nx, Ny, Nz = N

    span(ext) = ext === nothing ? one(FT) : (FT(ext[2]) - FT(ext[1]))
    Δx = span(x) / Nx
    Δy = span(y) / Ny
    Δz = span(z) / Nz

    return RectilinearGrid{TX, TY, TZ, FT}(Nx, Ny, Nz, FT(Δx), FT(Δy), FT(Δz))
end

# ---------------------------------------------------------------------------
# 3. The evaluation protocol
# ---------------------------------------------------------------------------

abstract type AbstractOperand end

# Numbers evaluate to themselves and are location-agnostic (never interpolated).
@inline evaluate(x::Number, i, j, k, grid) = x
location(::Number) = (Nothing, Nothing, Nothing)

# ---------------------------------------------------------------------------
# 4. Field (the leaf node)
# ---------------------------------------------------------------------------

struct Field{LX, LY, LZ, G <: RectilinearGrid, A <: AbstractArray} <: AbstractOperand
    grid::G
    data::A
end

location(::Field{LX, LY, LZ}) where {LX, LY, LZ} = (LX, LY, LZ)

Field{LX, LY, LZ}(grid::RectilinearGrid, data::AbstractArray) where {LX, LY, LZ} =
    Field{LX, LY, LZ, typeof(grid), typeof(data)}(grid, data)

"""
    Field{LX,LY,LZ}(grid)

Allocate a zeroed field at location `(LX, LY, LZ)`, sized for `grid`'s topology.
Use `Field{LX,LY,LZ}(grid, data)` to wrap an existing array instead.
"""
function Field{LX, LY, LZ}(grid::RectilinearGrid) where {LX, LY, LZ}
    TX, TY, TZ = topologies(grid)
    nx = npoints(TX, LX, grid.Nx)
    ny = npoints(TY, LY, grid.Ny)
    nz = npoints(TZ, LZ, grid.Nz)
    return Field{LX, LY, LZ}(grid, zeros(eltype(grid), nx, ny, nz))
end

# Topology-aware index wrapping: periodic wraps, bounded clamps (one-sided at the
# edge), flat collapses to 1. Because a flat axis always wraps to index 1, both
# δ (→ 0) and ℑ (→ identity) along it come out right with no special-casing.
@inline wrap_index(::Type{Periodic}, i, N) = mod1(i, N)
@inline wrap_index(::Type{Bounded},  i, N) = clamp(i, 1, N)
@inline wrap_index(::Type{Flat},     i, N) = 1

@inline function evaluate(f::Field, i, j, k, grid)
    TX, TY, TZ = topologies(grid)
    ii = wrap_index(TX, i, size(f.data, 1))
    jj = wrap_index(TY, j, size(f.data, 2))
    kk = wrap_index(TZ, k, size(f.data, 3))
    return @inbounds f.data[ii, jj, kk]
end

# ---------------------------------------------------------------------------
# 5. Index-level difference and interpolation kernels
#
# Dispatch on the *operand's* location along the worked axis; the result lands at
# the flipped location. Forward stencil moves Face → Center, backward Center → Face.
# ---------------------------------------------------------------------------

# --- x ---
@inline δx(::Type{Face},   i, j, k, grid, op) = evaluate(op, i+1, j, k, grid) - evaluate(op, i,   j, k, grid)
@inline δx(::Type{Center}, i, j, k, grid, op) = evaluate(op, i,   j, k, grid) - evaluate(op, i-1, j, k, grid)
@inline ℑx(::Type{Face},   i, j, k, grid, op) = (evaluate(op, i,   j, k, grid) + evaluate(op, i+1, j, k, grid)) / 2
@inline ℑx(::Type{Center}, i, j, k, grid, op) = (evaluate(op, i-1, j, k, grid) + evaluate(op, i,   j, k, grid)) / 2

# --- y ---
@inline δy(::Type{Face},   i, j, k, grid, op) = evaluate(op, i, j+1, k, grid) - evaluate(op, i, j,   k, grid)
@inline δy(::Type{Center}, i, j, k, grid, op) = evaluate(op, i, j,   k, grid) - evaluate(op, i, j-1, k, grid)
@inline ℑy(::Type{Face},   i, j, k, grid, op) = (evaluate(op, i, j,   k, grid) + evaluate(op, i, j+1, k, grid)) / 2
@inline ℑy(::Type{Center}, i, j, k, grid, op) = (evaluate(op, i, j-1, k, grid) + evaluate(op, i, j,   k, grid)) / 2

# --- z ---
@inline δz(::Type{Face},   i, j, k, grid, op) = evaluate(op, i, j, k+1, grid) - evaluate(op, i, j, k,   grid)
@inline δz(::Type{Center}, i, j, k, grid, op) = evaluate(op, i, j, k,   grid) - evaluate(op, i, j, k-1, grid)
@inline ℑz(::Type{Face},   i, j, k, grid, op) = (evaluate(op, i, j, k,   grid) + evaluate(op, i, j, k+1, grid)) / 2
@inline ℑz(::Type{Center}, i, j, k, grid, op) = (evaluate(op, i, j, k-1, grid) + evaluate(op, i, j, k,   grid)) / 2

# ---------------------------------------------------------------------------
# 6. Interpolation as a node (per-dimension location reconciliation)
# ---------------------------------------------------------------------------

# `shift?(from, to, ...)`: identity when locations match, else the right average.
@inline shiftx(::Type{L}, ::Type{L}, i, j, k, grid, op) where {L} = evaluate(op, i, j, k, grid)
@inline shiftx(::Type{Face},   ::Type{Center}, i, j, k, grid, op) = ℑx(Face,   i, j, k, grid, op)
@inline shiftx(::Type{Center}, ::Type{Face},   i, j, k, grid, op) = ℑx(Center, i, j, k, grid, op)

@inline shifty(::Type{L}, ::Type{L}, i, j, k, grid, op) where {L} = evaluate(op, i, j, k, grid)
@inline shifty(::Type{Face},   ::Type{Center}, i, j, k, grid, op) = ℑy(Face,   i, j, k, grid, op)
@inline shifty(::Type{Center}, ::Type{Face},   i, j, k, grid, op) = ℑy(Center, i, j, k, grid, op)

@inline shiftz(::Type{L}, ::Type{L}, i, j, k, grid, op) where {L} = evaluate(op, i, j, k, grid)
@inline shiftz(::Type{Face},   ::Type{Center}, i, j, k, grid, op) = ℑz(Face,   i, j, k, grid, op)
@inline shiftz(::Type{Center}, ::Type{Face},   i, j, k, grid, op) = ℑz(Center, i, j, k, grid, op)

struct InterpX{LX, LY, LZ, SRC, O} <: AbstractOperand; op::O; end
struct InterpY{LX, LY, LZ, SRC, O} <: AbstractOperand; op::O; end
struct InterpZ{LX, LY, LZ, SRC, O} <: AbstractOperand; op::O; end

location(::InterpX{LX, LY, LZ}) where {LX, LY, LZ} = (LX, LY, LZ)
location(::InterpY{LX, LY, LZ}) where {LX, LY, LZ} = (LX, LY, LZ)
location(::InterpZ{LX, LY, LZ}) where {LX, LY, LZ} = (LX, LY, LZ)

@inline evaluate(n::InterpX{LX, LY, LZ, SRC}, i, j, k, grid) where {LX, LY, LZ, SRC} =
    shiftx(SRC, LX, i, j, k, grid, n.op)
@inline evaluate(n::InterpY{LX, LY, LZ, SRC}, i, j, k, grid) where {LX, LY, LZ, SRC} =
    shifty(SRC, LY, i, j, k, grid, n.op)
@inline evaluate(n::InterpZ{LX, LY, LZ, SRC}, i, j, k, grid) where {LX, LY, LZ, SRC} =
    shiftz(SRC, LZ, i, j, k, grid, n.op)

"""
    interpolate(op, target::Tuple) -> AbstractOperand

Return a lazy node that evaluates `op` at location `target = (LX, LY, LZ)`, inserting a
half-cell average on each axis whose location differs from `target` and leaving matching
axes untouched (a `Number` is returned unchanged). This is the staggering primitive: e.g.
`compute(interpolate(H, (Face, Center, Center)))` moves an `aa`-node field onto the `acx`
velocity points. Materialize with [`compute`](@ref) / [`compute!`](@ref); to stagger a
whole expression instead of a single field, use [`@at`](@ref).
"""
# Bring `op` to `target = (LX, LY, LZ)` by nesting single-dim shifts (skipping any
# axis already at target). Nesting x-then-y gives the correct 4-point average,
# because the outer average evaluates the inner node at shifted indices.
interpolate(x::Number, target::Tuple) = x

function interpolate(op::AbstractOperand, target::Tuple)
    op = _interp_x(op, target[1])
    op = _interp_y(op, target[2])
    op = _interp_z(op, target[3])
    return op
end

function _interp_x(op, target)
    l = location(op)
    return l[1] === target ? op : InterpX{target, l[2], l[3], l[1], typeof(op)}(op)
end
function _interp_y(op, target)
    l = location(op)
    return l[2] === target ? op : InterpY{l[1], target, l[3], l[2], typeof(op)}(op)
end
function _interp_z(op, target)
    l = location(op)
    return l[3] === target ? op : InterpZ{l[1], l[2], target, l[3], typeof(op)}(op)
end

# ---------------------------------------------------------------------------
# 7. Operations: derivatives, binary/unary maps, and their front-ends
# ---------------------------------------------------------------------------

struct Derivative{LX, LY, LZ, DIR, O} <: AbstractOperand; op::O; end
location(::Derivative{LX, LY, LZ}) where {LX, LY, LZ} = (LX, LY, LZ)

@inline evaluate(d::Derivative{LX, LY, LZ, :x}, i, j, k, grid) where {LX, LY, LZ} =
    δx(flip(LX), i, j, k, grid, d.op) / grid.Δx
@inline evaluate(d::Derivative{LX, LY, LZ, :y}, i, j, k, grid) where {LX, LY, LZ} =
    δy(flip(LY), i, j, k, grid, d.op) / grid.Δy
@inline evaluate(d::Derivative{LX, LY, LZ, :z}, i, j, k, grid) where {LX, LY, LZ} =
    δz(flip(LZ), i, j, k, grid, d.op) / grid.Δz

# `∂?(L, op)`: produce a derivative whose result sits at `L`. To make the δ stencil
# land on the result location, the operand must sit at `flip` of that axis, so we
# interpolate it there first.
function ∂x(L::Tuple, op)
    Lx, Ly, Lz = L
    op2 = interpolate(op, (flip(Lx), Ly, Lz))
    return Derivative{Lx, Ly, Lz, :x, typeof(op2)}(op2)
end
function ∂y(L::Tuple, op)
    Lx, Ly, Lz = L
    op2 = interpolate(op, (Lx, flip(Ly), Lz))
    return Derivative{Lx, Ly, Lz, :y, typeof(op2)}(op2)
end
function ∂z(L::Tuple, op)
    Lx, Ly, Lz = L
    op2 = interpolate(op, (Lx, Ly, flip(Lz)))
    return Derivative{Lx, Ly, Lz, :z, typeof(op2)}(op2)
end

# Default (no `@at`): natural location = flip the operand's location on that axis.
∂x(op::AbstractOperand) = (l = location(op); ∂x((flip(l[1]), l[2], l[3]), op))
∂y(op::AbstractOperand) = (l = location(op); ∂y((l[1], flip(l[2]), l[3]), op))
∂z(op::AbstractOperand) = (l = location(op); ∂z((l[1], l[2], flip(l[3])), op))

# Binary / unary maps: interpolate every operand to L, then apply the scalar fn.
struct Operation{LX, LY, LZ, F, A} <: AbstractOperand
    f::F
    args::A
end
location(::Operation{LX, LY, LZ}) where {LX, LY, LZ} = (LX, LY, LZ)

@inline evaluate(o::Operation, i, j, k, grid) =
    o.f(map(a -> evaluate(a, i, j, k, grid), o.args)...)

function _binary(L::Tuple, f, args...)
    ops = map(a -> interpolate(a, L), args)
    return Operation{L[1], L[2], L[3], typeof(f), typeof(ops)}(f, ops)
end

# `a^n`: capture the exponent in a callable, then it's just a unary map.
_power(L::Tuple, op, n) = _binary(L, Base.Fix2(^, n), op)

# Front-ends without `@at` inherit a sensible default location (first operand's).
for op in (:+, :-, :*, :/)
    @eval Base.$op(a::AbstractOperand, b::AbstractOperand) = _binary(location(a), $op, a, b)
    @eval Base.$op(a::AbstractOperand, b::Number)          = _binary(location(a), $op, a, b)
    @eval Base.$op(a::Number, b::AbstractOperand)          = _binary(location(b), $op, a, b)
end
Base.:^(a::AbstractOperand, n::Number) = _power(location(a), a, n)
Base.sqrt(a::AbstractOperand)          = _binary(location(a), sqrt, a)

# ---------------------------------------------------------------------------
# 8. The `@at` macro: force every operation in `expr` to evaluate at `loc`
# ---------------------------------------------------------------------------

"""
    @at (LX, LY, LZ) expr

Build a lazy operation that evaluates `expr` at location `(LX, LY, LZ)`, inserting
interpolation wherever an operand's location doesn't match. Materialize with
`compute(op)` (or `compute!(field, op)`).
"""
macro at(loc, expr)
    return esc(insert_location(expr, loc))
end

# Walk the AST and thread `loc` into every recognized operator call. Function
# objects (`∂x`, `_binary`, ...) are embedded directly so the result needs no
# special scoping in the caller.
function insert_location(ex, loc)
    if ex isa Expr && ex.head === :call
        fsym = ex.args[1]
        rec  = Any[insert_location(a, loc) for a in ex.args[2:end]]
        if fsym === :∂x
            return Expr(:call, ∂x, loc, rec...)
        elseif fsym === :∂y
            return Expr(:call, ∂y, loc, rec...)
        elseif fsym === :∂z
            return Expr(:call, ∂z, loc, rec...)
        elseif fsym in (:+, :-, :*, :/)
            return Expr(:call, _binary, loc, getfield(Base, fsym), rec...)
        elseif fsym === :^
            return Expr(:call, _power, loc, rec...)
        elseif fsym === :sqrt
            return Expr(:call, _binary, loc, sqrt, rec...)
        else
            return Expr(:call, fsym, rec...)   # unknown call: recurse, no location
        end
    elseif ex isa Expr
        return Expr(ex.head, Any[insert_location(a, loc) for a in ex.args]...)
    else
        return ex
    end
end

# ---------------------------------------------------------------------------
# 9. Materialization
# ---------------------------------------------------------------------------

find_grid(f::Field)        = f.grid
find_grid(d::Derivative)   = find_grid(d.op)
find_grid(n::InterpX)      = find_grid(n.op)
find_grid(n::InterpY)      = find_grid(n.op)
find_grid(n::InterpZ)      = find_grid(n.op)
find_grid(o::Operation)    = _find_grid(o.args...)
_find_grid(a, rest...)     = a isa AbstractOperand ? find_grid(a) : _find_grid(rest...)

"""
    compute!(field, op)

Evaluate `op` into `field` (which must share `op`'s location). One scalar pass over
the array; the boundary is handled by the index wrapping in `evaluate(::Field, ...)`.
"""
function compute!(out::Field, op)
    grid = out.grid
    nx, ny, nz = size(out.data)
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        out.data[i, j, k] = evaluate(op, i, j, k, grid)
    end
    return out
end

"""
    compute(op) -> Field

Allocate a `Field` at `op`'s location and fill it by evaluating `op` everywhere.
"""
function compute(op::AbstractOperand)
    grid = find_grid(op)
    out  = Field{location(op)...}(grid)
    return compute!(out, op)
end

# ---------------------------------------------------------------------------
# 10. GPU support: make the operand tree `Adapt`-able
#
# Moving a node to a device adapts its leaf storage (the `Field` data array) and
# rebuilds the tree around it, preserving every location/type parameter. This is what
# lets `compute_ka!` push a node through a `CUDABackend` kernel: when the kernel is
# launched, the backend's adaptor converts each `CuArray` to a `CuDeviceArray` and
# recurses via these rules. `RectilinearGrid` (isbits) and `Number`s adapt to
# themselves through Adapt's identity fallback, so they need no rule.
# ---------------------------------------------------------------------------

Adapt.adapt_structure(to, f::Field{LX, LY, LZ}) where {LX, LY, LZ} =
    Field{LX, LY, LZ}(Adapt.adapt(to, f.grid), Adapt.adapt(to, f.data))

Adapt.adapt_structure(to, n::InterpX{LX, LY, LZ, SRC}) where {LX, LY, LZ, SRC} =
    (op = Adapt.adapt(to, n.op); InterpX{LX, LY, LZ, SRC, typeof(op)}(op))
Adapt.adapt_structure(to, n::InterpY{LX, LY, LZ, SRC}) where {LX, LY, LZ, SRC} =
    (op = Adapt.adapt(to, n.op); InterpY{LX, LY, LZ, SRC, typeof(op)}(op))
Adapt.adapt_structure(to, n::InterpZ{LX, LY, LZ, SRC}) where {LX, LY, LZ, SRC} =
    (op = Adapt.adapt(to, n.op); InterpZ{LX, LY, LZ, SRC, typeof(op)}(op))

Adapt.adapt_structure(to, d::Derivative{LX, LY, LZ, DIR}) where {LX, LY, LZ, DIR} =
    (op = Adapt.adapt(to, d.op); Derivative{LX, LY, LZ, DIR, typeof(op)}(op))

Adapt.adapt_structure(to, o::Operation{LX, LY, LZ}) where {LX, LY, LZ} =
    (args = map(a -> Adapt.adapt(to, a), o.args);
     Operation{LX, LY, LZ, typeof(o.f), typeof(args)}(o.f, args))

# ---------------------------------------------------------------------------
# TODO / natural extensions
#   * Stretched grids: store Δxᶜ / Δxᶠ coordinate vectors and pick spacing by the
#     result location, i.e. replace `grid.Δx` with `Δx(loc, grid, i)`.
#   * Halo regions + `fill_halo!` per topology, to replace index wrapping with
#     branch-free stencils (closer to Oceananigans, friendlier on GPU).
#   * Wrap `compute!` in a `KernelAbstractions.@kernel`; the `evaluate` tree is
#     type-stable and immutable, so it ports to GPU essentially unchanged.
# ---------------------------------------------------------------------------

end # module StaggeredGrids

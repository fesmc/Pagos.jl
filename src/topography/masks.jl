###############################################################
# Ice masks
###############################################################
#
# Three boolean cell-centred (`aa`) fields describe where a computation is meaningful:
#
#  - `is_ice`            — the cell holds ice (`H > H_min`).
#  - `is_ice_neighbour`  — the cell is ice-*free* but touches an ice cell: the one-cell
#                          ring the margin can advance into. Conceptually complementary to
#                          `TopographicMasks.is_margin` (the ring *inside* the ice) — but
#                          note `is_margin` is a pre-existing struct field that nothing in
#                          `src/` computes yet, so treat that pairing as a naming intent,
#                          not a working relationship, until `is_margin` gets a producer.
#  - `is_ice_allowed`    — a hard geometric constraint on where ice may exist at all.
#                          Not physical; imposed when a run must not grow ice outside a
#                          prescribed domain. Never derived from the state — it is set by
#                          the user and left alone by [`icemasks!`](@ref).
#
# Two distinct reasons to mask, and only the first is load-bearing:
#
#  1. **Containment.** Several constitutive quantities are not merely uninteresting on
#     ice-free cells, they are `NaN` — and the `NaN` *spreads into the ice*. The viscosity
#     is zero where there is no ice, `hlerp` averages reciprocals, so `τ_xy` at an ice-free
#     face is `NaN`, and the effective stress at `aa` interpolates that back one cell into
#     the ice, i.e. exactly onto the margin. Note this cannot be fixed by masking the
#     *strain rate*: `NaN * 0 == NaN`, so the kernel that touches the viscosity has to be
#     masked itself.
#  2. **Work avoidance**, which is the weaker motive here. These kernels are
#     bandwidth-bound and the mask is an extra load; on a GPU a divergent branch inside a
#     warp costs about the max of both paths, so the saving is real only where a whole warp
#     is uniformly ice-free. The construct that would genuinely save work at ice-sheet
#     coverage fractions is a compacted index list, not a branch. Masking is applied here
#     for correctness; do not expect a speedup from it without measuring.
#
# Chmy's own `FieldMask`/masked operators are deliberately *not* used: they cover only the
# derivative family (`left`/`right`/`δ`/`∂`/`∂²`/`∂k∂`/`divg`/`lapl`), and there is **no
# masked `lerp`/`hlerp`/`itp`** — so they cannot address the interpolation that produces the
# `NaN` above. They also carry a different semantics (multiply the operand by a mask weight
# before differencing, rather than skip), want eight float fields per mask on a 3D grid, and
# have no test coverage upstream in 0.1.26.

"""
$(TYPEDSIGNATURES)

Whether a computation is evaluated at a given node. Subtypes: [`NoMask`](@ref) (evaluate
everywhere) and [`IceMask`](@ref).

Masks are an **optional trailing argument** on the Chmy-native kernels, defaulting to
`NoMask()`, which dispatches away at compile time — so the unmasked path is branch-free and
identical to having no mask support at all.
"""
abstract type AbstractIceMask end

"""
$(TYPEDSIGNATURES)

Evaluate everywhere: the default. Compiles away entirely.
"""
struct NoMask <: AbstractIceMask end

"""
$(TYPEDSIGNATURES)

Evaluate only where the ice is (optionally, only inside a permitted domain).

A node is active when **any** cell it is built from is set in **any** of the `active`
fields, and — if `allowed` is given — **every** cell it is built from is permitted. The two
quantifiers differ on purpose: activity should extend to a face that has ice on one side
only (that face carries the flux that advances the margin), while a hard constraint must
block a face as soon as either side is out of bounds (so no mass crosses it, and none is
destroyed).

`active` and `allowed` are cell-centred (`aa`) boolean fields on the depth-integrated grid;
the node class the mask is queried at supplies the stencil of cells to reduce over, so one
mask serves every node class without per-location copies.

# Examples

```julia
IceMask(topo.mask.is_ice)                                       # strictly where ice is
IceMask(topo.mask.is_ice, topo.mask.is_ice_neighbour)           # ...plus the growth ring
IceMask(topo.mask.is_ice; allowed = topo.mask.is_ice_allowed)   # clipped to a fixed domain
```
"""
struct IceMask{F <: Tuple, A} <: AbstractIceMask
    active::F
    allowed::A
end
Adapt.@adapt_structure IceMask

# `active` is constrained to `<: Tuple` so that this varargs form, and not the struct's own
# two-positional-argument constructor, is what `IceMask(is_ice, is_ice_neighbour)` reaches.
# Without the constraint the second field silently becomes `allowed` instead of a second
# `active` field — which is a mask that quietly means something else entirely.
IceMask(active::AbstractArray...; allowed = nothing) = IceMask(active, allowed)

# The cells a node at horizontal location `(lx, ly)` is built from. This is Chmy's
# staggering convention — vertex `i` sits between centres `i - 1` and `i` — written out
# once, here, rather than re-derived at each call site. The vertical location is
# irrelevant: the masks are `grid2d` fields with no z dependence, so they are always read
# at `k = 1`.
@inline _mask_cells(::Center, ::Center, i, j) = ((i, j),)
@inline _mask_cells(::Vertex, ::Center, i, j) = ((i - 1, j), (i, j))
@inline _mask_cells(::Center, ::Vertex, i, j) = ((i, j - 1), (i, j))
@inline _mask_cells(::Vertex, ::Vertex, i, j) =
    ((i - 1, j - 1), (i, j - 1), (i - 1, j), (i, j))

# Is any of the `active` fields set at this one cell?
@inline _cell_set(fs::Tuple, c) = first(fs)[c..., 1] | _cell_set(Base.tail(fs), c)
@inline _cell_set(::Tuple{}, c) = false

@inline _any_cells_set(fs, cells::Tuple) = _cell_set(fs, first(cells)) |
                                           _any_cells_set(fs, Base.tail(cells))
@inline _any_cells_set(fs, ::Tuple{}) = false

@inline _all_cells_set(fs, cells::Tuple) = _cell_set(fs, first(cells)) &
                                           _all_cells_set(fs, Base.tail(cells))
@inline _all_cells_set(fs, ::Tuple{}) = true

# The recursion is a separate function from the `nothing` shortcut: sharing one name makes
# `(::Nothing, ::Tuple{})` match both the shortcut and the base case, with neither more
# specific — an ambiguity rather than a no-op.
@inline _all_allowed(::Nothing, cells) = true
@inline _all_allowed(f, cells) = _all_true(f, cells)

@inline _all_true(f, cells::Tuple) = f[first(cells)..., 1] & _all_true(f, Base.tail(cells))
@inline _all_true(f, ::Tuple{}) = true

"""
$(TYPEDSIGNATURES)

Whether the node at horizontal indices `(i, j)` and node class `loc` is active under `mask`:
**any** cell the node is built from is active, and (if `allowed` is set) every one of them is
permitted. `loc` is a full three-location tuple (one of the `NODE_*` constants); only its
horizontal part is consulted, since the masks have no vertical dependence.

This is the permissive rule, and it is the right one for anything whose value at a margin
node is physically meaningful: a mass flux (the face between ice and no ice is exactly the
one that advances the margin), a surface gradient, a velocity gradient, a strain rate. For
quantities that *interpolate* a field which vanishes off-ice — the viscosity — see
[`node_fully_active`](@ref).
"""
@inline node_active(::NoMask, loc, i, j) = true

@inline function node_active(m::IceMask, loc, i, j)
    cells = _mask_cells(loc[1], loc[2], i, j)
    return _any_cells_set(m.active, cells) & _all_allowed(m.allowed, cells)
end

"""
$(TYPEDSIGNATURES)

The strict counterpart of [`node_active`](@ref): **every** cell the node is built from must
be active (and permitted).

Needed wherever a node's value interpolates a quantity that is zero off-ice, the viscosity
being the case that matters: `hlerp` averages reciprocals, so a corner with one ice-free
neighbour gives `0/0 = NaN` under the permissive rule even though the harmonic mean of
`(η, 0)` is `0` in the limit — no stress is transmitted through a cell of zero viscosity.
The strict rule therefore delivers the physically correct value *and* avoids the `NaN`,
rather than trading one for the other. Callers do not choose between the two: the kernels
that need this apply it themselves to whatever mask they are given.
"""
@inline node_fully_active(::NoMask, loc, i, j) = true

@inline function node_fully_active(m::IceMask, loc, i, j)
    cells = _mask_cells(loc[1], loc[2], i, j)
    return _all_cells_set(m.active, cells) & _all_allowed(m.allowed, cells)
end

###############################################################
# Deriving the masks
###############################################################

@kernel inbounds = true function _is_ice!(is_ice, H, H_min, O)
    I = @index(Global, NTuple)
    I = I + O
    is_ice[I...] = H[I...] > H_min
end

@kernel inbounds = true function _is_ice_neighbour!(is_ice_neighbour, is_ice, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, k = I
    is_ice_neighbour[I...] = !is_ice[I...] &
                             (is_ice[i - 1, j, k] | is_ice[i + 1, j, k] |
                              is_ice[i, j - 1, k] | is_ice[i, j + 1, k])
end

"""
$(TYPEDSIGNATURES)

Recompute `mask.is_ice` and `mask.is_ice_neighbour` from the ice thickness, treating a cell
as ice-covered when `H > H_min`. `mask.is_ice_allowed` is **not** touched: it is a
user-imposed constraint, not a function of the state.

Runs two launches: `is_ice` is pointwise, `is_ice_neighbour` is a stencil over it and so
cannot share the pass. Only the interiors are meaningful — `is_ice_neighbour`'s halo ring
reads `is_ice` one cell further out than `is_ice`'s own sweep reached.

!!! note "Call this whenever the ice extent changes"
    The neighbour ring is one cell wide, which is exactly enough for one flux-divergence
    evaluation: in flux form mass can only enter a cell across a face carrying flux, so it
    can never appear more than one cell beyond the ice in a single evaluation, whatever the
    time step. A multi-stage time integrator evaluates the divergence several times per
    step, though, so it needs either a mask refresh between substages or a band as wide as
    the number of stages. (What the time step itself has to respect is positivity — not
    draining a cell below zero — which is the usual advective CFL condition, a separate
    constraint from the mask width.)
"""
function icemasks!(mask::TopographicMasks, H, rt::Runtime; H_min = 0)
    H_min = convert(eltype(H), H_min)
    rt.launch2d(rt.arch, rt.grid2d, _is_ice! => (mask.is_ice, H, H_min))
    rt.launch2d(rt.arch, rt.grid2d,
                _is_ice_neighbour! => (mask.is_ice_neighbour, mask.is_ice))
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level [`icemasks!`](@ref): derives the masks from `topo.thickness.ice`.
"""
icemasks!(topo::TopographicState, rt::Runtime; kwargs...) =
    icemasks!(topo.mask, topo.thickness.ice, rt; kwargs...)

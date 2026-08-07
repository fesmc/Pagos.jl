###############################################################
# Ice masks
###############################################################
#
# Three boolean cell-centred (`aa`) fields describe where a computation is meaningful:
#
#  - `is_ice`            — the cell holds ice (`H > H_min`).
#  - `is_ice_neighbour`  — ice-free but touches ice: the one-cell ring the margin can
#                          advance into.
#  - `is_ice_allowed`    — a hard, user-imposed geometric constraint on where ice may
#                          exist; never derived from the state, left alone by
#                          [`icemasks!`](@ref).
#
# @dev: `TopographyMasks.is_margin` is named to pair with `is_ice_neighbour` (the ring
# *inside* the ice), but nothing in `src/` computes it yet — treat the pairing as naming
# intent, not a working relationship, until it gets a producer.
#
# Masking matters for correctness, not just performance: several constitutive quantities
# are `NaN` off-ice (viscosity is zero there, and `hlerp` averages reciprocals), and the
# `NaN` spreads one cell into the ice through interpolation. Masking the strain rate instead
# doesn't help, since `NaN * 0 == NaN` — the kernel that touches viscosity has to be masked
# itself. Work avoidance is a secondary, weaker motive: these kernels are bandwidth-bound and
# a GPU warp only saves time when uniformly ice-free, so don't expect a speedup from masking
# alone without measuring.
#
# Chmy's own `FieldMask`/masked operators are not used here: they cover only the derivative
# family (no masked `lerp`/`hlerp`/`itp`, so they can't address the `NaN` above), multiply by
# a mask weight rather than skip, cost eight float fields per mask in 3D, and have no test
# coverage upstream (0.1.26).

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
struct IceMask{F<:Tuple,A} <: AbstractIceMask
    active::F
    allowed::A
end
Adapt.@adapt_structure IceMask

# `active` is constrained to `<: Tuple` so `IceMask(is_ice, is_ice_neighbour)` reaches this
# varargs constructor rather than the struct's own two-positional one — without it, the
# second field would silently become `allowed`, a mask that means something else entirely.
IceMask(active::AbstractArray...; allowed = nothing) = IceMask(active, allowed)

# The cells a node at horizontal location `(lx, ly)` is built from, per Chmy's staggering
# convention (vertex `i` sits between centres `i - 1` and `i`). Vertical location doesn't
# matter: the masks are `grid2d` fields with no z dependence, always read at `k = 1`.
@inline _mask_cells(::Center, ::Center, i, j) = ((i, j),)
@inline _mask_cells(::Vertex, ::Center, i, j) = ((i - 1, j), (i, j))
@inline _mask_cells(::Center, ::Vertex, i, j) = ((i, j - 1), (i, j))
@inline _mask_cells(::Vertex, ::Vertex, i, j) =
    ((i - 1, j - 1), (i, j - 1), (i - 1, j), (i, j))

# Is any of the `active` fields set at this one cell?
@inline _cell_set(fs::Tuple, c) = first(fs)[c..., 1] | _cell_set(Base.tail(fs), c)
@inline _cell_set(::Tuple{}, c) = false

@inline _any_cells_set(fs, cells::Tuple) =
    _cell_set(fs, first(cells)) | _any_cells_set(fs, Base.tail(cells))
@inline _any_cells_set(fs, ::Tuple{}) = false

@inline _all_cells_set(fs, cells::Tuple) =
    _cell_set(fs, first(cells)) & _all_cells_set(fs, Base.tail(cells))
@inline _all_cells_set(fs, ::Tuple{}) = true

# Kept separate from the `nothing` shortcut: sharing one name would make
# `(::Nothing, ::Tuple{})` match both, an ambiguity rather than a no-op.
@inline _all_allowed(::Nothing, cells) = true
@inline _all_allowed(f, cells) = _all_true(f, cells)

@inline _all_true(f, cells::Tuple) =
    f[first(cells)..., 1] & _all_true(f, Base.tail(cells))
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
    is_ice_neighbour[I...] =
        !is_ice[I...] &
        (is_ice[i-1, j, k] | is_ice[i+1, j, k] | is_ice[i, j-1, k] | is_ice[i, j+1, k])
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
function icemasks!(mask::TopographyMasks, H, rt::Runtime; H_min = 0)
    H_min = convert(eltype(H), H_min)
    rt.launch2d(rt.arch, rt.grid2d, _is_ice! => (mask.is_ice, H, H_min))
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _is_ice_neighbour! => (mask.is_ice_neighbour, mask.is_ice),
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level [`icemasks!`](@ref): derives the masks from `topo.thickness.ice`.
"""
icemasks!(topo::TopographicState, rt::Runtime; kwargs...) =
    icemasks!(topo.mask, topo.thickness.ice, rt; kwargs...)

###############################################################
# Where the momentum balance is well-posed
###############################################################
#
# Detached floating ice (an iceberg) has no basal drag and no membrane connection to the
# rest of the sheet, so nothing balances its driving stress: the SSA/DIVA system restricted
# to it is singular in its rigid-translation modes, and no solver setting fixes that — the
# problem posed there has no answer to converge to. This earns a dedicated mask rather than
# a docstring caveat for three reasons:
#
#  1. **It poisons global norms.** The stopping criterion is a max-norm, so a handful of
#     runaway faces can set `err` for the whole domain — on the AIS example this pinned the
#     scaled residual at 3.8e-3 from just 0.02% of the ice.
#  2. **It would poison the autotuner.** Duretz Eq. 21 estimates λ_min from a Rayleigh
#     quotient; a rigid-translation null mode drives λ_min → 0 and hence the damping
#     `c_damp` → 0 for the entire domain (`roadmaps/PT-autotune.md`, Phase 1.5).
#  3. **It costs nothing to act on.** Masks are already a per-call argument: momentum takes
#     `IceMask(is_momentum_solved)` while advection keeps `IceMask(is_ice, is_ice_neighbour)`
#     — bergs still advect and calve, they just stop being asked to satisfy a force balance.
#
# Connectivity is 4-way because that's what this discretization implies, not a convention: a
# corner-only contact does reach `N_xy` (at `ab`), but `_membrane_stress_staggered!` gates
# that term on `node_fully_active`, so a diagonally touching berg carries zero stress here.

# See [`momentum_mask!`](@ref) for why the in-place update below is race-safe.
@kernel inbounds = true function _grow_momentum_mask!(m, is_ice, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, k = I
    if is_ice[I...] & !m[I...]
        m[I...] = m[i-1, j, k] | m[i+1, j, k] | m[i, j-1, k] | m[i, j+1, k]
    end
end

@kernel inbounds = true function _seed_momentum_mask!(m, is_ice, is_grounded, O)
    I = @index(Global, NTuple)
    I = I + O
    m[I...] = is_ice[I...] & is_grounded[I...]
end

"""
$(TYPEDSIGNATURES)

Fill `is_momentum_solved`: the cells where the momentum balance is well-posed, i.e. ice
that is connected to grounded ice through ice. Detached icebergs are excluded — see the
source note above for why they are not merely uninteresting but unsolvable, and why
4-connectivity is the right rule for our stencil.

Pass the result to the momentum solver as its mask, in place of `is_ice`:

```julia
momentum_mask!(topo, rt)
pseudo_transient!(mech, cst, solver, rt, momentum, IceMask(topo.mask.is_momentum_solved))
```

and keep giving advection its own `IceMask(topo.mask.is_ice, topo.mask.is_ice_neighbour)`,
so detached ice still advects and calves — it is only excused from a force balance it
cannot satisfy. Nothing in the solver needs to know this mask exists; masks are already a
per-call argument.

`is_grounded` is the seed and is **read, not computed** (like `is_ice_allowed`, and unlike
`is_ice`): grounding is a flotation/grounding-line question that belongs to whatever
produced it, and this function works the same whether it came from a flotation criterion or
a sub-grid `f_grnd` threshold. Cells that are grounded but ice-free are ignored, so the
seed is `is_ice & is_grounded`.

# Implementation

Unlike every other mask here, connectivity is a *global* property and cannot be settled by
one local stencil pass. This is iterative label propagation — seed, then repeatedly grow
into ice neighbours until nothing changes — rather than a serial flood fill, so it stays one
`Launcher` sweep per iteration and runs unchanged on GPU. Growth is monotone (false → true
only), which makes the in-place update race-safe without a ping-pong buffer: a thread
reading a stale `false` just grows on the next sweep instead. The cost is `O(sweeps)`
boolean passes, paid once per change of ice extent (not per PT iteration, and not per time
step unless the extent moved): against the thousands of PT iterations of a single solve it
does not register.

Termination is two *consecutive* unchanged interior counts, not one. The `Launcher` sweeps
one halo ring beyond the interior, so the halo can be one sweep ahead of the interior;
stopping on the first unchanged count could therefore drop a cell that was about to be
reached from the ring. The second ghost ring is never written and so is never a source.

Throws if `maxsweeps` is exhausted rather than returning a partial mask: an under-grown mask
does not fail loudly, it silently freezes real, connected shelf ice at zero velocity — the
same class of quiet wrong answer that `is_momentum_solved` exists to eliminate.
"""
function momentum_mask!(
    is_momentum_solved,
    is_ice,
    is_grounded,
    rt::Runtime;
    maxsweeps = 4 * sum(size(rt.grid2d, Center())),
)
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _seed_momentum_mask! => (is_momentum_solved, is_ice, is_grounded),
    )

    m = asarray(is_momentum_solved)
    count_prev = sum(m)
    quiet = 0
    for _ = 1:maxsweeps
        rt.launch2d(
            rt.arch,
            rt.grid2d,
            _grow_momentum_mask! => (is_momentum_solved, is_ice),
        )
        count_now = sum(m)
        quiet = count_now == count_prev ? quiet + 1 : 0
        quiet == 2 && return nothing
        count_prev = count_now
    end
    throw(
        ErrorException(
            "momentum_mask!: connectivity did not converge in maxsweeps = $maxsweeps sweeps. " *
            "Raise `maxsweeps` — the count needed scales with the longest ice path across the " *
            "domain, which a serpentine geometry can make much longer than its diameter. " *
            "Returning the partial mask instead would silently freeze connected ice at zero " *
            "velocity, so this throws.",
        ),
    )
end

"""
$(TYPEDSIGNATURES)

State-level [`momentum_mask!`](@ref): writes `topo.mask.is_momentum_solved` from
`topo.mask.is_ice` and `topo.mask.is_grounded`. Call after [`icemasks!`](@ref) (it reads
`is_ice`) and after whatever sets `is_grounded`.
"""
momentum_mask!(topo::TopographicState, rt::Runtime; kwargs...) = momentum_mask!(
    topo.mask.is_momentum_solved,
    topo.mask.is_ice,
    topo.mask.is_grounded,
    rt;
    kwargs...,
)

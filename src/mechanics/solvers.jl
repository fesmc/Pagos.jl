###############################################################
# Solvers
##############################################################

"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the dynamics solver via [`velocity`](@ref).

# Available subtypes:
"""
abstract type AbstractMomentumSolver end

# TODO
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via energy minimization.
"""
struct OptimMomentumSolver <: AbstractMomentumSolver
end

# TODO
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via an iterative linear solver (e.g., CG, GMRES).
"""
struct IterativeMomentumSolver <: AbstractMomentumSolver
end

"""
$(TYPEDSIGNATURES)

TODO: Solve the ice dynamics via wavelet methods.
"""
struct WaveletMomentumSolver <: AbstractMomentumSolver
end

# TODO
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via convolutional neural network (IGM style).
"""
struct ConvolutionalMomentumSolver <: AbstractMomentumSolver
end

# TODO
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via a transient solver (e.g., explicit time-stepping).
"""
struct TransientMomentumSolver <: AbstractMomentumSolver
end

"""
$(TYPEDSIGNATURES)

An abstract type to multiple-dispatch whether and how [`PseudoTransientSolver`](@ref)
updates `material.viscosity_depthaveraged` from the current velocity iterate during the
PT loop, following the codebase's "dispatch, not `if`/`else`" convention
(`roadmaps/pagos.md`). [`NoViscosityContinuation`](@ref) (the default) makes the solver
treat viscosity as a fixed input, exactly as every solver behaved before this type existed.
[`GlenViscosityContinuation`](@ref) enables Sandip et al. (2024) Eq. 8's log-space
continuation.

!!! note "Standalone, not integrated with `AbstractFlowLaw`/`AbstractCreep`"
    `GlenViscosityContinuation` implements Glen's law directly from the strain-rate
    invariant (Sandip Eq. 3), independent of `src/material/flow_law.jl`/`creep.jl`. Those
    are built around a *stress*-driven creep formulation (`creep(σ_e, law)`), which would
    need an inner per-cell nonlinear solve to use here (`σ_e = 2ηε̇_e` depends on the very
    `η` being solved for) — not what Sandip's direct, closed-form update does. This is
    flagged, not resolved, in `roadmaps/PT-autotune.md`: the codebase now has two
    parallel notions of "Glen's law" (this one and the stress-driven one), plus two
    *dead* strain-rate-driven copies (`src/legacy/flow_law.jl`,
    `src/thermodynamics/viscosity.jl`, neither `include`d) that were the closer relatives
    of this one and went unused instead of being reused. Making this coherent — one flow-
    law abstraction, or an explicit, documented reason for two — is future work.
"""
abstract type AbstractViscosityContinuation end

"""
$(TYPEDSIGNATURES)

No viscosity continuation: `material.viscosity_depthaveraged` is a fixed input for the
whole PT solve, untouched by the solver. The default for [`PseudoTransientSolver`](@ref),
and the only behaviour prior to [`AbstractViscosityContinuation`](@ref) existing.
"""
struct NoViscosityContinuation <: AbstractViscosityContinuation end

"""
$(TYPEDSIGNATURES)

Glen-law viscosity continuation (Sandip et al. 2024, Eq. 3 and 8): at every PT iteration,
`material.viscosity_depthaveraged` is updated from the current effective strain rate
`ε̇_e` (Sandip Eq. 3, `η = A^(-1/n)·(ε̇_e² + ε̇₀²)^((1-n)/(2n))/2`, regularized by the
strain-rate floor `strainrate_reg = ε̇₀`) and relaxed toward its own previous value in
log-space (Eq. 8, `η_new = exp(θ_μ·log(η_raw) + (1-θ_μ)·log(η_old))`) rather than applied
directly — direct application is what Sandip found could diverge the PT iteration before
the nonlinear viscosity has had time to relax, `θ_μ` controlling how much time.

# Fields:
 - `n_glen`: Glen's flow-law exponent (default `3`, matching [`GlenNyeCreep`](@ref)'s).
 - `theta_mu`: log-space relaxation weight, `0 < theta_mu ≤ 1` (`1` disables relaxation:
   `η_new` is `η_raw` outright, no memory of `η_old`).
 - `strainrate_reg`: strain-rate floor `ε̇₀` preventing the `ε̇_e = 0` singularity
   (`η → ∞` as `ε̇_e → 0` in Glen's law); must be chosen relative to the problem's own
   strain-rate scale, there is no unit-independent default that means the same thing
   everywhere — the constructor requires it rather than guessing.

Reads `material.rate_factor_depthaveraged` as the (prescribed, not thermally coupled) rate
factor `A` — see that field's docstring.
"""
struct GlenViscosityContinuation{T<:AbstractFloat} <: AbstractViscosityContinuation
    n_glen::T
    theta_mu::T
    strainrate_reg::T
end

function GlenViscosityContinuation(T::Type{<:AbstractFloat} = Float64;
    n_glen = 3, theta_mu = 0.1, strainrate_reg,
)
    return GlenViscosityContinuation{T}(T(n_glen), T(theta_mu), T(strainrate_reg))
end

"""
$(TYPEDSIGNATURES)

An abstract type to multiple-dispatch how [`PseudoTransientSolver`](@ref) chooses its
pseudo-time step `Δτ`, following the same "dispatch, not `if`/`else`" convention as
[`AbstractViscosityContinuation`](@ref). [`ViscosityPseudoTimeStep`](@ref) (the default) is
Sandip et al. (2024) Eq. 7, `Δτ ∝ 1/η`; [`GershgorinPseudoTimeStep`](@ref) bounds the
spectral radius of the *actual* SSA/DIVA residual operator instead.

Only read by the Chmy-native, C-grid staggered [`pseudo_transient!`](@ref); the collocated
solver always uses its own global scalar [`pseudo_dt`](@ref).
"""
abstract type AbstractPseudoTimeStep end

"""
$(TYPEDSIGNATURES)

Sandip et al. (2024) Eq. 7: `Δτ = dtau_scaling · ρ dx dy / (4 (1 + muB) ndim2 · η_face)`,
evaluated per grid point from the local depth-averaged viscosity ([`pseudo_dt!`](@ref)).
The default, and the only behaviour prior to [`AbstractPseudoTimeStep`](@ref) existing.

!!! warning "Ignores basal drag and the ice mask"
    The formula bounds the membrane-stress (diffusive) part of the operator only. It knows
    nothing about the basal-drag term `β u`, which is the *dominant* diagonal contribution
    under grounded ice, and it interpolates the viscosity across the ice margin, so an
    ice-free cell's placeholder viscosity throttles `Δτ` on the margin faces. Both are why
    a real-geometry solve needs [`GershgorinPseudoTimeStep`](@ref) — see its docstring and
    `roadmaps/PT-autotune.md`.
"""
struct ViscosityPseudoTimeStep <: AbstractPseudoTimeStep end

"""
$(TYPEDSIGNATURES)

Pseudo-time step from a Gershgorin bound on the spectral radius of the SSA/DIVA residual
operator (Duretz et al. 2026, Eq. 20 — `roadmaps/PT-autotune.md` Phase 2): the explicit
stability limit is `Δτ ≤ 2/λ_max`, and `λ_max` is bounded by the largest absolute row sum
of the (mass-scaled) operator, which for the `u`-equation at an `acx` face is

```
Λ_x = [ 8(P₋+P₊)/dx² + 4(P₋+P₊)/(dx dy) + 2(Q₋+Q₊)/dy² + 2(Q₋+Q₊)/(dx dy) + β_face ] / (ρ H_face)
```

with `P = ηH` at the two adjacent `aa` cells, `Q = hlerp(η)·lerp(H)` at the two adjacent
`ab` corners (exactly the coefficients the membrane-stress [`strainrate!`](@ref) builds),
`β_face = lerp(β)` and `H_face = lerp(H)`. The `y`-equation is the mirror image.
`Δτ = cfl · 2 / Λ`.

Three things this buys over [`ViscosityPseudoTimeStep`](@ref):

 1. **The basal-drag term is in the bound.** `β/(ρH)` dominates `λ_max` under grounded ice
    (`β_eff` up to ~5e13 Pa s m⁻¹ on real geometry), and a `Δτ` that ignores it diverges
    outright as soon as friction is solved for rather than prescribed.
 2. **Correct magnitude for the membrane part.** With uniform `η`, `H`, no drag and
    `dx = dy`, this gives `Δτ = ρ dx²/(16 η)` — which is what Sandip's Eq. 7 reduces to at
    `muB = 0`, `ndim2 = 4.1`. The `PseudoTransientSolver` default `muB = 1e2` therefore
    shrinks the step by ~101×, i.e. costs ~101× the iterations, for no stability benefit on
    a depth-integrated balance (`muB` is a *bulk*-viscosity ratio: it belongs to the
    compressible/full-Stokes pressure step of Räss et al. 2020, not to SSA/DIVA, which has
    no pressure unknown).
 3. **Mask-consistent at the margin.** Every coefficient is evaluated through the same
    [`AbstractIceMask`](@ref) the membrane stress uses, so an ice-free neighbour
    contributes `0` (it transmits no stress) instead of throttling `Δτ` with whatever
    placeholder viscosity that cell happens to hold.

# Fields:
 - `cfl`: safety factor on `2/λ_max`, `0 < cfl ≤ 1` (default `0.9`). The Gershgorin bound
   is an upper bound on `λ_max`, so `cfl = 1` is already conservative in exact arithmetic;
   the margin is there for the nonlinear case (`GlenViscosityContinuation` moves `η`
   between the `Δτ` evaluation and its use).

!!! note "Still computed once per solve"
    Like [`ViscosityPseudoTimeStep`](@ref), this is evaluated once before the PT loop, from
    the viscosity and friction fields as they stand then. Re-estimating it *during* the
    loop (needed when `η` evolves under viscosity continuation) is the re-estimation cadence
    item of `roadmaps/PT-autotune.md` Phase 2, not yet implemented.
"""
struct GershgorinPseudoTimeStep{T<:AbstractFloat} <: AbstractPseudoTimeStep
    cfl::T
end

GershgorinPseudoTimeStep(T::Type{<:AbstractFloat} = Float64; cfl = 0.9) =
    GershgorinPseudoTimeStep{T}(T(cfl))

"""
$(TYPEDSIGNATURES)

An abstract type to multiple-dispatch what [`pseudo_transient!`](@ref) compares against
`solver.abstol`, following the same "dispatch, not `if`/`else`" convention as
[`AbstractViscosityContinuation`](@ref). [`VelocityIncrement`](@ref) (the default) is the
max-norm velocity change per iteration; [`ScaledResidual`](@ref) is the momentum residual
normalized by the driving-stress scale.

Only read by the Chmy-native, C-grid staggered [`pseudo_transient!`](@ref).
"""
abstract type AbstractPTConvergence end

"""
$(TYPEDSIGNATURES)

Stop on the max-norm velocity increment `max|u_new - u_old|` (units of `u`). The default,
and the only behaviour prior to [`AbstractPTConvergence`](@ref) existing.

!!! warning "Silently reports convergence on a stiff sub-domain"
    The increment is `Δτ · r(u)`, so it goes to zero wherever `Δτ` is small — converged or
    not. On real geometry that is not academic: an ice shelf has no basal drag, so it
    relaxes purely diffusively and its per-iteration increment is orders of magnitude below
    the grounded ice's from the very first iteration. A tolerance that the grounded ice has
    to work to reach is one the shelf satisfies while still at ~0 velocity, and the solve
    returns `converged = true` with the shelves empty. Use [`ScaledResidual`](@ref) when the
    domain mixes drag regimes.
"""
struct VelocityIncrement <: AbstractPTConvergence end

"""
$(TYPEDSIGNATURES)

Stop on the momentum residual, nondimensionalized by the driving-stress scale:

```
err = max|r(u)| / max|τ_d / (ρ H_face)|
```

where `r(u) = (∇·N - τ_b - τ_d)/(ρH)` is what [`dotvel!`](@ref) already writes into
`solver.residual_x`/`residual_y`, and the normalization is the same quantity evaluated with
the membrane and basal terms dropped — i.e. the rate the driving stress alone would produce.
`err` is therefore dimensionless and `O(1)` at a zero initial guess, so `abstol` reads as
"fraction of the driving-stress forcing left unbalanced" (`1e-3`–`1e-6` are sensible), not
as a velocity.

!!! warning "`abstol` changes units when you select this"
    With [`VelocityIncrement`](@ref) `abstol` is in m s⁻¹; here it is dimensionless. This is
    exactly the semantic change `roadmaps/PT-autotune.md` Phase 1 deferred rather than force
    on every existing test — hence a dispatch type rather than a change of meaning in place.

Normalized by the driving stress rather than by the *initial* residual on purpose: a
transient run re-solves from the previous time step's velocity, where the initial residual
is already small, and a relative-reduction criterion would then silently demand many orders
more accuracy than the first solve got.
"""
struct ScaledResidual <: AbstractPTConvergence end

"""
$(TYPEDSIGNATURES)

An abstract type to multiple-dispatch whether [`PseudoTransientSolver`](@ref) recomputes
`stress.base_x`/`base_y` from the current velocity iterate during the PT loop, following
the same "dispatch, not `if`/`else`" convention as [`AbstractViscosityContinuation`](@ref).
[`ActiveFrictionUpdate`](@ref) (the default) is the ordinary basal friction law, recomputed
every iteration. [`NoFrictionUpdate`](@ref) leaves `stress.base_x`/`base_y` untouched for
the whole solve, so whatever was written there beforehand (e.g. a prescribed basal stress
field) stays fixed instead of being overwritten from `beta_eff * velocity` — i.e. the
friction law is bypassed.
"""
abstract type AbstractFrictionUpdate end

"""
$(TYPEDSIGNATURES)

The ordinary basal friction law: `stress.base_x`/`base_y` are recomputed every PT
iteration from `friction.beta_eff * velocity.base_{x,y}` via [`basalstress!`](@ref). The
default for [`PseudoTransientSolver`](@ref), and the only behaviour prior to
[`AbstractFrictionUpdate`](@ref) existing.
"""
struct ActiveFrictionUpdate <: AbstractFrictionUpdate end

"""
$(TYPEDSIGNATURES)

Bypass the basal friction law: `stress.base_x`/`base_y` are never touched by the PT loop,
so they stay at whatever value they were set to before the solve — a fixed, prescribed
basal stress rather than one derived from a friction law and the current velocity.
"""
struct NoFrictionUpdate <: AbstractFrictionUpdate end

"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via a pseudo-transient (PT) solver following Sandip et al. (2024).
The velocity field is relaxed in pseudo-time until the momentum balance is satisfied,
which only requires local (stencil) operations. All work arrays live on the backend of
the grid the solver is constructed from, so the same code runs on CPU and GPU.

Construct from a grid, overriding parameters selectively via keyword arguments:

```julia
solver = PseudoTransientSolver(grid)
solver = PseudoTransientSolver(grid; maxiter = 500, abstol = 1e-10)
```

# Fields:
 - `ndim1`, `ndim2`, `ndim3`: numerical-dimensionality constants of the PT time step
   (1D, 2D, 3D stencils). Read only by [`ViscosityPseudoTimeStep`](@ref).
 - `muB`: bulk-to-shear viscosity ratio entering the PT time step. Read only by
   [`ViscosityPseudoTimeStep`](@ref) — see [`GershgorinPseudoTimeStep`](@ref) for why the
   default `1e2` is ~101× too conservative on a depth-integrated balance.
 - `theta_v`: relaxation weight of the velocity update.
 - `gamma`: damping coefficient of the pseudo-transient rate (Sandip et al. 2024,
   Eq. 12–14; Frankel 1950). The rate is accumulated as
   `dv_new = (1 - gamma) * dv_old + r(u)`, so `gamma = 1` discards all memory and
   recovers the plain (undamped, first-order) PT iteration — the default, and the only
   value exercised by the collocated solver. `0 < gamma < 1` turns the iteration into a
   damped wave, which is what buys sub-quadratic iteration-count scaling; picking a
   good value currently needs the same manual tuning as `theta_v`. Automatically
   selecting it from spectral estimates is the Duretz et al. (2026) autotuning work
   (`roadmaps/PT-autotune.md`), not yet implemented.
 - `abstol`: convergence tolerance on whatever `convergence` measures — a velocity change
   per iteration (m s⁻¹) for [`VelocityIncrement`](@ref), a dimensionless residual for
   [`ScaledResidual`](@ref).
 - `maxiter`: maximum number of PT iterations.
 - `ncheck`: check convergence every `ncheck` iterations. The check is the only
   operation that forces the host to wait for the device (see
   [Asynchronous kernel launches](@ref)), so on GPU a larger value (10–50) keeps
   the launch pipeline busy at the cost of up to `ncheck - 1` extra iterations.
 - `printout_every`: print a convergence monitor every `printout_every` iterations
   (silent by default).
 - `dtau_scaling`: safety scaling of the PT time step.
 - `velocity_x_old`, `velocity_y_old`: previous velocity iterate.
 - `velocity_x_dt`, `velocity_y_dt`: pseudo-transient velocity rate (the damping
   accumulator described under `gamma`).
 - `dtau_x`, `dtau_y`: pseudo-time step, held at the same node class as the velocity
   component it advances. Only the Chmy-native, C-grid staggered
   [`pseudo_transient!`](@ref) fills these per grid point (from the local viscosity, via
   [`pseudo_dt!`](@ref)); the collocated solver keeps its single global scalar (from
   [`pseudo_dt`](@ref)) and never reads these fields, which is why they may look unused
   in that code path.
 - `residual_x`, `residual_y`: the raw (undamped) pseudo-transient rate `r(u)/(ρH)`,
   written unconditionally by [`dotvel!`](@ref) regardless of `gamma` — a genuinely
   separate buffer from `velocity_x_dt`/`velocity_y_dt` because damping's whole point is
   for the accumulator to *not* equal the raw rate, so the raw rate needs its own home to
   remain readable. The Chmy-native [`pseudo_transient!`](@ref) reduces these into the
   `residual` diagnostic it returns; unused by the collocated solver, like `dtau_x`/`dtau_y`.
 - `pseudo_timestep`: an [`AbstractPseudoTimeStep`](@ref), default
   [`ViscosityPseudoTimeStep`](@ref). Selects how `dtau_x`/`dtau_y` are filled. Only read
   by the Chmy-native, C-grid staggered [`pseudo_transient!`](@ref).
 - `convergence`: an [`AbstractPTConvergence`](@ref), default [`VelocityIncrement`](@ref).
   Selects what `abstol` is compared against — **and therefore what `abstol` means**. Only
   read by the Chmy-native, C-grid staggered [`pseudo_transient!`](@ref).
 - `viscosity_continuation`: an [`AbstractViscosityContinuation`](@ref), default
   [`NoViscosityContinuation`](@ref). Only read by the Chmy-native, C-grid staggered
   [`pseudo_transient!`](@ref); the collocated solver never updates viscosity regardless
   of this field's value.
 - `friction_update`: an [`AbstractFrictionUpdate`](@ref), default
   [`ActiveFrictionUpdate`](@ref). Only read by the Chmy-native, C-grid staggered
   [`pseudo_transient!`](@ref); the collocated solver always recomputes basal stress
   inline regardless of this field's value. Set to [`NoFrictionUpdate`](@ref) to bypass
   the basal friction law and hold `stress.base_x`/`base_y` fixed at whatever was written
   there before the solve (e.g. a prescribed basal stress field).

!!! note "Two type parameters for the work arrays, not one"
    On a [`StaggeredGrid`](@ref) the x- and y-velocity work arrays are Chmy `Field`s at
    different locations (`acx`/`acy`), which are different concrete types (location is
    part of a `Field`'s type). A single `M` shared by all four fields — as a plain-array
    `RegularGrid` state can get away with, since `Matrix{T}` carries no location — would
    reject that combination outright, hence the `MX`/`MY` split below. They collapse to
    the same type on a `RegularGrid`, so that path is unaffected.
"""
struct PseudoTransientSolver{T<:AbstractFloat, MX, MY, PT<:AbstractPseudoTimeStep, CV<:AbstractPTConvergence, VC<:AbstractViscosityContinuation, FU<:AbstractFrictionUpdate} <: AbstractMomentumSolver
    ndim1::T
    ndim2::T
    ndim3::T
    muB::T
    theta_v::T
    gamma::T
    abstol::T
    maxiter::Int
    ncheck::Int
    printout_every::Int
    dtau_scaling::T
    velocity_x_old::MX
    velocity_y_old::MY
    velocity_x_dt::MX
    velocity_y_dt::MY
    dtau_x::MX
    dtau_y::MY
    residual_x::MX
    residual_y::MY
    pseudo_timestep::PT
    convergence::CV
    viscosity_continuation::VC
    friction_update::FU
end
Adapt.@adapt_structure PseudoTransientSolver

function PseudoTransientSolver(grid::RegularGrid;
    ndim1 = 2.1,
    ndim2 = 4.1,
    ndim3 = 6.1,
    muB = 1e2,
    theta_v = 0.6,
    gamma = 1,
    abstol = 1e-8,
    maxiter = 100,
    ncheck = 1,
    printout_every = typemax(Int),
    dtau_scaling = 1,
    pseudo_timestep::AbstractPseudoTimeStep = ViscosityPseudoTimeStep(),
    convergence::AbstractPTConvergence = VelocityIncrement(),
    viscosity_continuation::AbstractViscosityContinuation = NoViscosityContinuation(),
    friction_update::AbstractFrictionUpdate = ActiveFrictionUpdate(),
)
    T = eltype(grid.x)
    backend = get_backend(grid.x)
    w() = KernelAbstractions.zeros(backend, T, grid.nx, grid.ny)
    return PseudoTransientSolver(
        T(ndim1), T(ndim2), T(ndim3), T(muB),
        T(theta_v), T(gamma), T(abstol), Int(maxiter), Int(ncheck),
        Int(printout_every), T(dtau_scaling), w(), w(), w(), w(), w(), w(), w(), w(),
        pseudo_timestep, convergence, viscosity_continuation, friction_update,
    )
end

"""
$(TYPEDSIGNATURES)

Build a [`PseudoTransientSolver`](@ref) of Chmy `Field`s for the Chmy-native, C-grid
staggered [`pseudo_transient!`](@ref). The velocity-iterate/rate work arrays
(`velocity_x_old`, `velocity_y_old`, `velocity_x_dt`, `velocity_y_dt`) are placed at
`acx`/`acy` on `grid.grid2d`, matching the velocity components they mirror
(`mech.velocity.x`/`y`).

Requires a depth-averaged grid (`grid.grid2d === grid.grid`, i.e. `nz == 1`), checked at
construction rather than left to fail at first solve: DIVA's vertical-shear integral is
Phase 3 future work (`roadmaps/chmy.md`), so today's Chmy-native `pseudo_transient!` only
supports the depth-averaged case, exactly as the collocated method does in practice
(it reads `view(velocity.x, :, :, 1)`).
"""
function PseudoTransientSolver(grid::StaggeredGrid;
    ndim1 = 2.1,
    ndim2 = 4.1,
    ndim3 = 6.1,
    muB = 1e2,
    theta_v = 0.6,
    gamma = 1,
    abstol = 1e-8,
    maxiter = 100,
    ncheck = 1,
    printout_every = typemax(Int),
    dtau_scaling = 1,
    halo = 1,
    pseudo_timestep::AbstractPseudoTimeStep = ViscosityPseudoTimeStep(),
    convergence::AbstractPTConvergence = VelocityIncrement(),
    viscosity_continuation::AbstractViscosityContinuation = NoViscosityContinuation(),
    friction_update::AbstractFrictionUpdate = ActiveFrictionUpdate(),
)
    grid.grid2d === grid.grid || throw(ArgumentError(
        "PseudoTransientSolver(::StaggeredGrid) requires a depth-averaged grid " *
        "(grid.grid2d === grid.grid, i.e. nz == 1); DIVA's vertical shear integral is " *
        "not yet ported (roadmaps/chmy.md, Phase 3). Build the grid without a " *
        "`layering` argument, e.g. `StaggeredGrid(T, lx, ly, dx, dy)`."))

    (; arch) = grid
    g = grid.grid2d
    T = eltype(g)
    acx() = _field(arch, g, NODE_ACX, T, halo)
    acy() = _field(arch, g, NODE_ACY, T, halo)
    return PseudoTransientSolver(
        T(ndim1), T(ndim2), T(ndim3), T(muB),
        T(theta_v), T(gamma), T(abstol), Int(maxiter), Int(ncheck),
        Int(printout_every), T(dtau_scaling),
        acx(), acy(), acx(), acy(), acx(), acy(), acx(), acy(),
        pseudo_timestep, convergence, viscosity_continuation, friction_update,
    )
end

"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via a direct linear solver (e.g., sparse LU factorization).

# Improvements over `LegacyLinearMomentumSolver2D`:
 1. dynamics is a type parameter → dispatch on loop1!/loop2! is fully static, no need to thread a runtime `dynamics` argument through populate_vectors! layers.
 2. SparseMatrixCSC pre-allocated at construction (fixed sparsity pattern). The hot path writes directly to A.nzval via a precomputed COO→nzval index map, eliminating the sparse(Ai, Aj, Av) allocation on every solve.
 3. Single AI type parameter (i_idx and j_idx are always the same kind).
 4. VT/MT/PI type parameters for vector/matrix/perm arrays so the struct can hold GPU arrays (CuVector, CuSparseMatrix) without code changes. The populate_vectors! kernels are written with KernelAbstractions and run on whichever backend owns lsd.u.
"""
struct LinearMomentumSolver2D{
    DYN <: AbstractMomentumBalance,
    T   <: AbstractFloat,
    VT  <: AbstractVector,        # float vector type (u, u0, b)
    MT,                            # sparse matrix type (SparseMatrixCSC or CuSparseMatrix)
    PI  <: AbstractVector{Int},   # perm index vector type
    AI,
} <: AbstractMomentumSolver
    dynamics::DYN
    nx::Int
    ny::Int
    dxdx_::T
    dydy_::T
    dxdy_::T
    u::VT
    u0::VT
    b::VT
    A::MT
    perm::PI                      # COO fill order → A.nzval index
    i_idx::AI
    j_idx::AI
    solver_cache::Ref{Any}        # holds a cached direct solver (Nothing or backend-specific)
end

###############################################################
# Dispatch functions
###############################################################

# Functions to calculate velocity
function calc_F_integral(visc_eff, H_ice, f_ice, zeta_aa, n)
    # TODO: not yet implemented.
    error("calc_F_integral is not yet implemented")
end

# function velocity( solver::LinearSolver, dynamics::DIVA)
# end

function vertically_integrated_viscosity!(N, H, μ)
    N .= H .* μ
    return nothing
end

function stagger()
end

function fill_pattern!(Ai, Aj, nx, ny, i_idx, j_idx)
    k = 0
    for i in 1:nx, j in 1:ny
        im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
        nr = _ij2n_ux(i, j, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(ip1, j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(i,   j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(im1, j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(i,   jp1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(i,   jm1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(i,   j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(ip1, j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(ip1, jm1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(i,   jm1, nx, ny)
    end
    for i in 1:nx, j in 1:ny
        im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
        nr = _ij2n_uy(i, j, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(i,   jp1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(i,   j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(i,   jm1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(ip1, j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(im1, j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(i,   jp1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(i,   j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(im1, jp1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(im1, j,   nx, ny)
    end
    return nothing
end

function coo_to_nzval_idx(Ai, Aj, A::SparseMatrixCSC)
    perm = Vector{Int}(undef, length(Ai))
    for k in eachindex(Ai)
        c  = Aj[k]
        r  = Ai[k]
        lo = A.colptr[c]
        hi = A.colptr[c+1] - 1
        perm[k] = searchsortedfirst(A.rowval, r, lo, hi, Base.Order.Forward)
    end
    return perm
end

function LinearMomentumSolver2D(grid::RegularGrid, dynamics::DYN, backend = CPU()) where {DYN <: AbstractMomentumBalance}
    T      = eltype(grid.x)
    nx, ny = grid.nx, grid.ny
    dx, dy = T(grid.dx), T(grid.dy)
    dxdx_  = 1 / (dx * dx)
    dydy_  = 1 / (dy * dy)
    dxdy_  = 1 / (dx * dy)
    n_sprs = 18 * nx * ny  # 9 nonzeros/row × 2 equations (ux, uy)
    n_u    = 2 * nx * ny

    i_idx = PeriodicIndexing(1, nx)
    j_idx = PeriodicIndexing(1, ny)

    # Pattern and permutation are always computed on CPU (one-time cost).
    Ai_cpu   = zeros(Int, n_sprs)
    Aj_cpu   = zeros(Int, n_sprs)
    fill_pattern!(Ai_cpu, Aj_cpu, nx, ny, i_idx, j_idx)
    A_cpu    = sparse(Ai_cpu, Aj_cpu, ones(T, n_sprs), n_u, n_u)
    perm_cpu = coo_to_nzval_idx(Ai_cpu, Aj_cpu, A_cpu)

    # Allocate live arrays on the target backend.
    # For GPU (e.g. CUDABackend()), KernelAbstractions.zeros returns CuVector and
    # the sparse matrix should be adapted via CUDA.CUSPARSE.CuSparseMatrixCSC(A_cpu).
    u    = KernelAbstractions.zeros(backend, T,   n_u)
    u0   = KernelAbstractions.zeros(backend, T,   n_u)
    b    = KernelAbstractions.zeros(backend, T,   n_u)
    perm = KernelAbstractions.zeros(backend, Int, n_sprs)
    perm .= perm_cpu   # works for both CPU (no-op copy) and GPU (H→D transfer)

    # A stays as SparseMatrixCSC for the CPU default; for GPU, adapt before passing
    # to the inner constructor, e.g.:
    #   A = CUDA.CUSPARSE.CuSparseMatrixCSC(A_cpu)
    return LinearMomentumSolver2D(dynamics, nx, ny, dxdx_, dydy_, dxdy_, u, u0, b, A_cpu, perm, i_idx, j_idx, Ref{Any}(nothing))
end

@kernel function _assemble_ux!(nzval, perm, u0, b, nx, ny,
                               N, N_ab, ux, taud_acx, β_acx, β_acy,
                               dxdx_, dydy_, dxdy_, i_idx, j_idx, dynamics)
    i, j = @index(Global, NTuple)
    im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
    nr = _ij2n_ux(i, j, nx, ny)
    @inbounds u0[nr] = ux[i, j]
    @inbounds b[nr]  = taud_acx[i, j]
    v = loop1_coeffs(im1, i, ip1, jm1, j, jp1, dxdx_, dydy_, dxdy_,
                     N, N_ab, β_acx, β_acy, dynamics)
    k0 = 9 * ((i - 1) * ny + (j - 1))
    @inbounds for s in 1:9
        nzval[perm[k0 + s]] = v[s]
    end
end

@kernel function _assemble_uy!(nzval, perm, u0, b, nx, ny,
                               N, N_ab, uy, taud_acy, β_acx, β_acy,
                               dxdx_, dydy_, dxdy_, i_idx, j_idx, dynamics)
    i, j = @index(Global, NTuple)
    im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
    nr = _ij2n_uy(i, j, nx, ny)
    @inbounds u0[nr] = uy[i, j]
    @inbounds b[nr]  = taud_acy[i, j]
    v = loop2_coeffs(im1, i, ip1, jm1, j, jp1, dxdx_, dydy_, dxdy_,
                     N, N_ab, β_acx, β_acy, dynamics)
    k0 = 9 * nx * ny + 9 * ((i - 1) * ny + (j - 1))
    @inbounds for s in 1:9
        nzval[perm[k0 + s]] = v[s]
    end
end

"""
$(TYPEDSIGNATURES)

Assemble the linear system of the [`LinearMomentumSolver2D`](@ref) from the current dynamic
state `dyn_now`: fill the sparse matrix `A` (viscosity and basal-drag coefficients) and the
right-hand side `b` (driving stress), and seed the initial guess `u0` from the current
velocity. Mutates the solver buffers in place; the low-level method takes the unpacked fields
directly.
"""
function populate_vectors!(lsd::LinearMomentumSolver2D, dyn_now)
    (; N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy) = dyn_now
    populate_vectors!(lsd, N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy)
    return nothing
end

function populate_vectors!(
    lsd::LinearMomentumSolver2D,
    N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy,
)
    (; A, perm, u0, b, i_idx, j_idx, nx, ny, dxdx_, dydy_, dxdy_, dynamics) = lsd
    backend  = get_backend(lsd.u)
    kern_ux  = _assemble_ux!(backend)
    kern_uy  = _assemble_uy!(backend)

    nzval = nonzeros(A)
    kern_ux(
        nzval, perm, u0, b, nx, ny,
        N, N_ab, ux, taud_acx, β_acx, β_acy,
        dxdx_, dydy_, dxdy_, i_idx, j_idx, dynamics;
        ndrange = (nx, ny),
    )
    kern_uy(
        nzval, perm, u0, b, nx, ny,
        N, N_ab, uy, taud_acy, β_acx, β_acy,
        dxdx_, dydy_, dxdy_, i_idx, j_idx, dynamics;
        ndrange = (nx, ny),
    )
    return nothing
end

function LinearSolve.LinearProblem(lsd::LinearMomentumSolver2D)
    return LinearProblem(lsd.A, lsd.b; u0 = lsd.u)
end

function velocity!(lsd::LinearMomentumSolver2D)
    if lsd.solver_cache[] === nothing
        F = lu(lsd.A)
        lsd.solver_cache[] = F
        ldiv!(lsd.u, F, lsd.b)
    else
        _velocity_cached!(lsd.solver_cache[], lsd)
    end
    return nothing
end

# Function barrier: typed on F so lu! and ldiv! dispatch statically.
function _velocity_cached!(F, lsd::LinearMomentumSolver2D)
    lu!(F, lsd.A)
    ldiv!(lsd.u, F, lsd.b)
    return nothing
end

function velocity!(ux, uy, lsd::LinearMomentumSolver2D)
    u = lsd.u
    (; nx, ny) = lsd
    @inbounds for i in 1:nx, j in 1:ny
        ux[i, j] = u[_ij2n_ux(i, j, nx, ny)]
        uy[i, j] = u[_ij2n_uy(i, j, nx, ny)]
    end
    return nothing
end
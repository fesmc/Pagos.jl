###############################################################
# Solvers
##############################################################

"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the dynamics solver via [`velocity`](@ref).

# Available subtypes:
"""
abstract type AbstractMomentumSolver end

# @dev TODO: not yet implemented.
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via energy minimization.
"""
struct OptimMomentumSolver <: AbstractMomentumSolver end

# @dev TODO: not yet implemented.
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via an iterative linear solver (e.g., CG, GMRES).
"""
struct IterativeMomentumSolver <: AbstractMomentumSolver end

# @dev TODO: not yet implemented.
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via wavelet methods.
"""
struct WaveletMomentumSolver <: AbstractMomentumSolver end

# @dev TODO: not yet implemented.
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via convolutional neural network (IGM style).
"""
struct ConvolutionalMomentumSolver <: AbstractMomentumSolver end

# @dev TODO: not yet implemented.
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via a transient solver (e.g., explicit time-stepping).
"""
struct TransientMomentumSolver <: AbstractMomentumSolver end

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
 - `strainrate_reg`: strain-rate floor `ε̇₀` (yr⁻¹) preventing the `ε̇_e = 0` singularity
   (`η → ∞` as `ε̇_e → 0` in Glen's law); must be chosen relative to the problem's own
   strain-rate scale, there is no unit-independent default that means the same thing
   everywhere — the constructor requires it rather than guessing.

   Being an absolute-scale floor rather than a ratio, this is one of the few numbers that
   a change of time unit silently reinterprets: a value picked when the code ran in
   seconds regularizes ~3.16e7 times more weakly read as yr⁻¹. Ice deforms at order
   `1e-4`–`1e-2` yr⁻¹, so a floor well below that leaves Glen's law in charge.

Reads `material.rate_factor_depthaveraged` as the (prescribed, not thermally coupled) rate
factor `A` — see that field's docstring.
"""
struct GlenViscosityContinuation{T<:AbstractFloat} <: AbstractViscosityContinuation
    n_glen::T
    theta_mu::T
    strainrate_reg::T
end

function GlenViscosityContinuation(
    T::Type{<:AbstractFloat} = Float64;
    n_glen = 3,
    theta_mu = 0.1,
    strainrate_reg,
)
    return GlenViscosityContinuation{T}(T(n_glen), T(theta_mu), T(strainrate_reg))
end

"""
$(TYPEDSIGNATURES)

DIVA's Glen-law viscosity continuation: the same flow law as
[`GlenViscosityContinuation`](@ref), applied **per layer** on the column grid.

The difference is which effective strain rate feeds it. `GlenViscosityContinuation` uses
Eq. (12), the depth-averaged invariant, and writes the single field
`material.viscosity_depthaveraged`. This one uses Eq. (13) — the same invariant *plus* the
vertical-shear terms `¼(u_z² + v_z²)`, with `u_z` diagnosed from the basal stress via
Eq. (21) — and so produces a genuinely depth-varying `µ(z)` in `material.viscosity`. The
depth-averaged field is then derived from it by [`depthaverage!`](@ref), rather than
computed independently: `µ̄` must be the average of the `µ(z)` the shear was built from, or
the membrane and shear terms describe different ice.

Reads `material.rate_factor` (the column rate factor `A(z)`), not
`rate_factor_depthaveraged`.

Fields are as [`GlenViscosityContinuation`](@ref)'s: `n_glen`, `theta_mu`,
`strainrate_reg`.

!!! note "This is the only writer of `µ` under DIVA"
    Making it an [`AbstractViscosityContinuation`](@ref) rather than a separate mechanism
    keeps one dispatch point for "who owns the viscosity" (`roadmaps/chmy.md`, Phase 3,
    decision 12). Pairing `DIVAMomentumBalance` with `GlenViscosityContinuation` instead
    would leave `µ(z)` untouched and the shear terms stale — which is why
    [`diva_update!`](@ref) drives this directly rather than relying on the solver's
    continuation slot alone.

!!! note "`NoViscosityContinuation` leaves `µ(z)`/`µ̄` mutually unchecked, by design"
    Under [`NoViscosityContinuation`](@ref), `material.viscosity` and
    `material.viscosity_depthaveraged` are both prescribed inputs — like
    `rate_factor`/`rate_factor_depthaveraged` already are — and nothing here verifies that
    `µ̄` is actually the depth average of `µ(z)`. That consistency is the responsibility of
    whatever populates `MechanicState`'s material fields once topography, dynamics,
    thermodynamics and material are wired together as one model; it is not this solver's
    job to guess at it or guard against it. Exactly the same stance as decision 7's
    `β_eff = 0` contract: the caller owns the input, the docstring states the contract, no
    runtime check.
"""
struct DIVAViscosityContinuation{T<:AbstractFloat} <: AbstractViscosityContinuation
    n_glen::T
    theta_mu::T
    strainrate_reg::T
end

function DIVAViscosityContinuation(
    T::Type{<:AbstractFloat} = Float64;
    n_glen = 3,
    theta_mu = 0.1,
    strainrate_reg,
)
    return DIVAViscosityContinuation{T}(T(n_glen), T(theta_mu), T(strainrate_reg))
end

"""
$(TYPEDSIGNATURES)

[`BlatterPattynMomentumBalance`](@ref)'s viscosity continuation: the same per-layer Glen law
as [`DIVAViscosityContinuation`](@ref), writing `material.viscosity` (`µ(z)`) on `rt.grid`
from BP's effective strain rate ([`effective_strainrate_bp!`](@ref), Eq. 3), but from the
**actual** `velocity.x_dz`/`y_dz` rather than DIVA's Eq. 21 diagnosis from `τ_b` — BP already
carries a real 3D velocity, so there is nothing to diagnose.

Unlike [`DIVAViscosityContinuation`](@ref), `material.viscosity_depthaveraged` is *not*
derived here: nothing on the BP path reads it (`roadmaps/blatter-pattyn.md`, §1.3), so it is
left as whatever degenerate allocation the state constructor gave it, mirroring
`viscosity_integral_1`/`_2` and `rate_factor_depthaveraged` staying untouched on this path.

Reads `material.rate_factor` (the column rate factor `A(z)`), not
`rate_factor_depthaveraged` — the same choice as `DIVAViscosityContinuation`.

Fields are as [`GlenViscosityContinuation`](@ref)'s: `n_glen`, `theta_mu`, `strainrate_reg`.
"""
struct BPViscosityContinuation{T<:AbstractFloat} <: AbstractViscosityContinuation
    n_glen::T
    theta_mu::T
    strainrate_reg::T
end

function BPViscosityContinuation(
    T::Type{<:AbstractFloat} = Float64;
    n_glen = 3,
    theta_mu = 0.1,
    strainrate_reg,
)
    return BPViscosityContinuation{T}(T(n_glen), T(theta_mu), T(strainrate_reg))
end

###############################################################
# DIVA depth-integrated-viscosity update
###############################################################

"""
$(TYPEDSIGNATURES)

How often the DIVA chain `µ(z) → F₁/F₂ → β_eff → µ̄` is re-evaluated *during* a
pseudo-transient solve.

The chain is a fixed point — `µ` depends on `u_z`, which depends on `τ_b`, which depends on
`β_eff`, which depends on `F₂`, which depends on `µ` — and Robinson et al. (2022) prescribe
only that the quantities come "from the previous iteration", not how often that iteration
should be. Hence a strategy rather than a hard-coded cadence.

**DIV** is *depth-integrated viscosity*, the paper's own decomposition of the DIVA acronym
and exactly what the chain recomputes. Deliberately not named `…ViscosityUpdate`, which
would sit one word away from the unrelated [`NoViscosityContinuation`](@ref) on the same
solver.

Subtypes: [`NoDIVUpdate`](@ref) (default), [`PeriodicDIVUpdate`](@ref).
"""
abstract type AbstractDIVUpdate end

"""
$(TYPEDSIGNATURES)

Never refresh the DIVA chain during the solve: `β_eff`, `µ(z)`, `µ̄` and `F₁`/`F₂` are held
at whatever the caller put there, and the nonlinearity is carried by the outer timestep.
The default.

!!! warning "The caller owns the chain, and a zero `β_eff` is frictionless sliding"
    Under `NoDIVUpdate` the solver never evaluates the chain — not even once before the
    loop. A caller that skips [`diva_update!`](@ref) gets `friction.beta_eff` at its
    allocation default of zero, i.e. **no basal drag at all**, with no error raised. This is
    a deliberate contract (`roadmaps/chmy.md`, Phase 3, decision 7): the solver does what it
    says and no redundant work, exactly like the halo-filling contract documented on
    [`pseudo_transient!`](@ref).

!!! note "What `converged` means here"
    With no in-loop refresh, a converged solve is a converged *frozen-`β_eff`* problem, not
    a converged DIVA fixed point. The DIVA nonlinearity is resolved across outer timesteps,
    not within one `pseudo_transient!` call.
"""
struct NoDIVUpdate <: AbstractDIVUpdate end

"""
$(TYPEDSIGNATURES)

Re-evaluate the DIVA chain every `n_update` pseudo-transient iterations (`n_update = 1`
refreshes every iteration, the most faithful reading of the fixed point and the most
expensive; larger values amortize the cost).

# Fields
 - `n_update`: refresh period in PT iterations, `≥ 1`.

!!! note "A refresh also re-derives `Δτ`"
    `β_eff` is an input to the Gershgorin bound behind `solver.dtau_x`/`dtau_y`, computed
    once before the loop. A `β_eff` that grows under a stale bound makes that bound
    optimistic and the iteration can diverge, so every refresh re-runs
    [`pseudo_dt!`](@ref) (`roadmaps/chmy.md`, Phase 3, decision 8). This preserves the
    invariant `pseudo_dt!` already states for the mask and `friction_update`: the bound must
    describe the operator actually being iterated.
"""
struct PeriodicDIVUpdate <: AbstractDIVUpdate
    n_update::Int
    function PeriodicDIVUpdate(n_update::Integer = 1)
        n_update ≥ 1 || throw(
            ArgumentError(
                "PeriodicDIVUpdate requires n_update ≥ 1, got $n_update. Use NoDIVUpdate() to " *
                "disable in-loop refreshes entirely.",
            ),
        )
        return new(Int(n_update))
    end
end

"""
$(TYPEDSIGNATURES)

An abstract type to multiple-dispatch how [`PseudoTransientSolver`](@ref) chooses its
pseudo-time step `Δτ`, following the same "dispatch, not `if`/`else`" convention as
[`AbstractViscosityContinuation`](@ref). [`GershgorinPseudoTimeStep`](@ref) (the default)
bounds the spectral radius of the SSA/DIVA residual operator; [`ViscosityPseudoTimeStep`](@ref)
is Sandip et al. (2024) Eq. 7, `Δτ ∝ 1/η`.

Each subtype owns the constants its own formula needs, so no solver carries a parameter
only one `Δτ` rule reads.
"""
abstract type AbstractPseudoTimeStep end

"""
$(TYPEDSIGNATURES)

Sandip et al. (2024) Eq. 7: `Δτ = dtau_scaling · ρ dx dy / (4 (1 + muB) ndim · η_face)`,
evaluated per grid point from the local depth-averaged viscosity ([`pseudo_dt!`](@ref)).
The only behaviour prior to [`AbstractPseudoTimeStep`](@ref) existing.

# Fields:
 - `muB`: bulk-to-shear viscosity ratio (default `1e2`, Sandip's value).
 - `ndim`: numerical-dimensionality constant of the stencil (default `4.1`, the 2D value —
   this solver is depth-integrated, so the 1D and 3D constants Sandip also lists have
   nowhere to be used).

!!! warning "Ignores basal drag and the ice mask"
    The formula bounds the membrane-stress (diffusive) part of the operator only. It knows
    nothing about the basal-drag term `β u`, which is the *dominant* diagonal contribution
    under grounded ice, and it interpolates the viscosity across the ice margin, so an
    ice-free cell's placeholder viscosity throttles `Δτ` on the margin faces. Both are why
    a real-geometry solve needs [`GershgorinPseudoTimeStep`](@ref) — the default — and why
    this is kept mainly as the published formula to measure against.
"""
struct ViscosityPseudoTimeStep{T<:AbstractFloat} <: AbstractPseudoTimeStep
    muB::T
    ndim::T
end

ViscosityPseudoTimeStep(T::Type{<:AbstractFloat} = Float64; muB = 1e2, ndim = 4.1) =
    ViscosityPseudoTimeStep{T}(T(muB), T(ndim))

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
    With [`VelocityIncrement`](@ref) `abstol` is in m yr⁻¹; here it is dimensionless. This is
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

An abstract type to multiple-dispatch how [`PseudoTransientSolver`](@ref) obtains the two
parameters of the pseudo-transient iteration — the pseudo-time step `Δτ` and the damping
`γ` — following the same "dispatch, not `if`/`else`" convention as
[`AbstractViscosityContinuation`](@ref). [`FixedTuning`](@ref) (the default) carries them
as hand-set numbers; [`AutotunedDynamicRelaxation`](@ref) derives both from spectral
estimates during the solve (Duretz et al. 2026).

"""
abstract type AbstractPTTuning end

"""
$(TYPEDSIGNATURES)

Hand-set iteration parameters, held fixed for the whole solve: `Δτ` comes from
`solver.pseudo_timestep` alone, scaled by `theta_v`, and the damping is `gamma`. The
default, and the only behaviour prior to [`AbstractPTTuning`](@ref) existing.

Also the baseline [`AutotunedDynamicRelaxation`](@ref) is measured against — `gamma = 1`
is the undamped iteration the speedup numbers in `roadmaps/PT-autotune.md` are quoted
relative to, which is why hand-setting survives at all.

# Fields:
 - `theta_v`: relaxation weight of the velocity update (default `0.6`).
 - `gamma`: damping coefficient of the pseudo-transient rate (Sandip et al. 2024,
   Eq. 12–14; Frankel 1950). The rate is accumulated as
   `dv_new = (1 - gamma) * dv_old + r(u)`, so `gamma = 1` (the default) discards all
   memory and recovers the plain, undamped, first-order iteration; `0 < gamma < 1` turns
   it into a damped wave, which is what buys sub-quadratic iteration-count scaling.
   Finding a good value is a manual scan — the thing
   [`AutotunedDynamicRelaxation`](@ref) exists to replace.
"""
struct FixedTuning{T<:AbstractFloat} <: AbstractPTTuning
    theta_v::T
    gamma::T
end

FixedTuning(T::Type{<:AbstractFloat} = Float64; theta_v = 0.6, gamma = 1) =
    FixedTuning{T}(T(theta_v), T(gamma))

"""
$(TYPEDSIGNATURES)

Dynamic relaxation with automatically tuned `Δτ` and damping (Duretz et al. 2026;
`roadmaps/PT-autotune.md` Phase 2, which carries the derivation and the measurements).
Both are derived from spectral estimates of the *preconditioned* momentum operator rather
than scanned by hand, and re-derived every `cadence` iterations so they follow a viscosity
that evolves during the solve.

Writing the mass-scaled residual [`dotvel!`](@ref) forms as `r̃(u) = b̃ - Ã u` and `Λ` for
the Gershgorin absolute row sum of `Ã` that [`GershgorinPseudoTimeStep`](@ref) already
computes per face, the preconditioner is `M = diag(Λ)`. That makes `λ_max(M⁻¹Ã) ≤ 1` hold
*exactly* (the preconditioned row sums are 1), so `λ_max` costs no reduction and only
`λ_min` is measured — by a Rayleigh quotient over consecutive iterates the loop already
holds, `λ_min ≈ |Δuᵀ Δr̃| / (Δuᵀ M Δu)` (Duretz Eq. 21). Hence the requirement, checked at
construction, that `pseudo_timestep` be a [`GershgorinPseudoTimeStep`](@ref).

Duretz Eq. 19/20 give `c = c_damp·2√λ_min` and `Δτ = c_CFL·2/√λ_max` separately, but the
damping the loop applies is `γ = c·Δτ` and the damped iteration's exact stability condition
is `Δτ²λ_max ≤ 2(2 - γ)` — a `Δτ` chosen as if `γ = 0` sits outside it. Solving the pair
jointly at `λ_max = 1`, with `d = 2·c_damp·√λ_min`:

```
Δτ = -c_CFL²·d + √(c_CFL⁴·d² + 4·c_CFL²),    γ = d·Δτ  (always < 2)
```

`c_CFL` is `pseudo_timestep.cfl`, not a second knob; `cfl = 0.99` or tighter is intended
here (Duretz use `≲ 0.999`).

The first `cadence` iterations are a warm-up: undamped at `Δτ = 1/λ_max`, since the DR step
is unstable without the damping it is paired with and there is no initial `λ_min` to derive
that from. Its amplification `1 - λ̂` annihilates the `λ_max` mode, leaving `Δu` dominated by
the `λ_min` mode the first quotient must see.

# Fields:
 - `c_damp`: safety factor on `c = 2√λ_min`, `[0.5, 1]` per Duretz (default `0.8`). A
   Rayleigh quotient over-estimates `λ_min`, so `c_damp = 1` errs toward over-damping;
   measured, `0.5` and `0.8` are within ~10% and both beat `1.0` and `0.3`.
 - `cadence`: re-estimate every `cadence` iterations (default `50`), each re-estimation also
   refilling `Δτ` from a fresh Gershgorin bound. Duretz use ~100; tighter here because the
   first estimate is also what ends the warm-up. One re-estimation costs two reductions and
   a kernel — the cost of one `ncheck` check, which defaults to every iteration.

`theta_v` and `gamma` have no existence outside [`FixedTuning`](@ref), so there is nothing
here to leave unread: this type simply does not carry them.
[`pseudo_transient!`](@ref) returns the final `damping` and `lambda_min`.
"""
struct AutotunedDynamicRelaxation{T<:AbstractFloat} <: AbstractPTTuning
    c_damp::T
    cadence::Int
end

AutotunedDynamicRelaxation(
    T::Type{<:AbstractFloat} = Float64;
    c_damp = 0.8,
    cadence = 50,
) = AutotunedDynamicRelaxation{T}(T(c_damp), Int(cadence))

"""
$(TYPEDSIGNATURES)

Reject a `tuning`/`pseudo_timestep` pair whose spectral normalization does not hold, at
construction rather than as plausible-looking garbage at solve time. Only
[`AutotunedDynamicRelaxation`](@ref) constrains its partner.
"""
_check_tuning(::AbstractPTTuning, ::AbstractPseudoTimeStep) = nothing

_check_tuning(::AutotunedDynamicRelaxation, ::GershgorinPseudoTimeStep) = nothing

_check_tuning(::AutotunedDynamicRelaxation, pt::AbstractPseudoTimeStep) = throw(
    ArgumentError(
        "AutotunedDynamicRelaxation requires pseudo_timestep = GershgorinPseudoTimeStep(...), " *
        "got $(typeof(pt)). The autotuner's spectral estimates are normalized by the " *
        "Gershgorin row sum (λ_max ≤ 1 by construction); no such normalization holds for " *
        "another Δτ rule, so the derived Δτ and damping would be wrong by an unknown factor.",
    ),
)

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

###############################################################
# Vertical treatment (Blatter-Pattyn only)
###############################################################

"""
$(TYPEDSIGNATURES)

An abstract type to multiple-dispatch how a [`MomentumBalance3D`](@ref) solve advances the
vertical-shear divergence `∂z(µ ∂z u)` — explicitly, like every other term, or implicitly
down each column — following the same "dispatch, not `if`/`else`" convention as
[`AbstractPseudoTimeStep`](@ref). [`ExplicitVertical`](@ref) (the default) is the whole of
`roadmaps/blatter-pattyn.md` Phase 1; [`ImplicitVertical`](@ref) is Phase 2.

Read **only** on the [`MomentumBalance3D`](@ref) path: SSA and DIVA have no vertical operator
to treat (DIVA integrates it out analytically through `F₁`/`F₂`), so the depth-averaged
[`PseudoTransientSolver`](@ref) constructor rejects anything but `ExplicitVertical()` rather
than accepting a setting it would silently ignore.
"""
abstract type AbstractVerticalTreatment end

"""
$(TYPEDSIGNATURES)

Advance `∂z(µ ∂z u)` explicitly, together with the membrane terms: the velocity update is the
plain [`pseudo_vel!`](@ref) step `u ← u + θ_v Δτ dv`, and `Δτ` is bounded by the **full**
Gershgorin row sum, vertical term and bed-drag term included. The default, and the only
behaviour prior to [`AbstractVerticalTreatment`](@ref) existing.

The cost is the aspect-ratio penalty of `roadmaps/blatter-pattyn.md` §2: `Λ_vert/Λ_horiz` has
a median of ~400 and a 99th percentile of ~2e5 on 8 km Antarctic geometry, so `Δτ` is set by
the vertical operator almost everywhere and by the thinnest columns in particular
(`Λ_vert ∝ 1/H²`). [`ImplicitVertical`](@ref) removes exactly that penalty; this stays as the
reference it is verified against.
"""
struct ExplicitVertical <: AbstractVerticalTreatment end

"""
$(TYPEDSIGNATURES)

Advance `∂z(µ ∂z u)` **implicitly** down each column — vertical line relaxation,
`roadmaps/blatter-pattyn.md` Phase 2 — while the membrane terms stay explicit. `Δτ` is then
bounded by the horizontal operator alone, i.e. by the same row sum SSA and DIVA use, and the
aspect-ratio penalty of [`ExplicitVertical`](@ref) disappears rather than being square-rooted
by the damping.

Read as a preconditioner, which is what makes it a drop-in: the residual
[`dotvel!`](@ref) writes is unchanged (still the true, fully explicit `r̃(u^k)`), and only the
*step* changes, from `Δu = θ_v Δτ dv` to

```
(diag(Λ_horiz/ρ̃) + Ã_v) Δu = θ_v · scale · dv
```

with `Ã_v` the exact discrete vertical operator — the tridiagonal built from the same `R±`,
`δ±`, `Δ_k` coefficients the Gershgorin bound uses, so the two cannot drift apart. Solved
per column by the Thomas algorithm, one work-item per horizontal face, embarrassingly
parallel over `(i, j)` and GPU-friendly. Non-uniform sigma spacing is coefficients only —
[`CorrectedVerticalLayering`](@ref) needs no special path.

Two wins, not one. The vertical diffusion leaves the explicit spectrum, and so does the
**bed drag**: `β` enters the bottom row of the tridiagonal exactly, rather than the `β/Δζ₁`
term the explicit bound has to carry. On 8 km Antarctic geometry 87 % of the residual left
after 600 explicit iterations sits at `k = 1`, so that second win is the larger one.

!!! warning "Verified on clean geometry; does not yet converge on real Antarctic geometry"
    Same fixed point as [`ExplicitVertical`](@ref) and an iteration count flat in `nz` (1100 →
    1420 over `nz ∈ {4…32}`, against 1720 → 25 380 explicit) — on uniform slabs and on
    synthetically masked, laterally varying cases. On the real 8 km Antarctic restart the
    solve instead *cycles*: down to `err ~ 3e-2`, a burst to `1e6`–`1e7`, recovery, repeat.
    The cause is open and tracked in `roadmaps/blatter-pattyn.md`, Phase 2 ("recurring
    bursts"), which also records what has already been ruled out (the `cfl`/`λ_max` margin,
    and the `λ_min = 1` clamp). Prefer [`ExplicitVertical`](@ref) on real geometry.

# Fields
 - `thomas_x`, `thomas_y`: the Thomas back-substitution coefficients `c'`, one column field
   per velocity component. The forward sweep's `b'` and `d'` never outlive one layer (`d'` is
   written into the velocity field itself), so these two arrays are the entire extra
   footprint — allocated here rather than on [`PseudoTransientSolver`](@ref) so that an
   [`ExplicitVertical`](@ref) solve pays nothing for a strategy it does not use.

Construct from the same [`StaggeredGrid`](@ref) the solver is built on:

```julia
solver = PseudoTransientSolver(grid, BlatterPattynMomentumBalance();
                               vertical_treatment = ImplicitVertical(grid))
```

!!! note "The autotuner stays valid, with the implicit factor folded into `M`"
    [`AutotunedDynamicRelaxation`](@ref) needs `λ_max(M⁻¹Ã) ≤ 1`, which under
    [`ExplicitVertical`](@ref) holds because `M = diag(Λ)` is the full Gershgorin diagonal.
    Here `M = diag(Λ_horiz/ρ̃) + Ã_v` — the preconditioner the line solve actually applies —
    and `M - Ã = diag(Λ_horiz/ρ̃) - Ã_h` is exactly the horizontal Gershgorin defect, so the
    bound survives verbatim. Folding it in costs nothing: `Δuᵀ M Δu = scale · Δuᵀ dv` by
    construction of the step, so the `M`-inner product is a reduction over two fields the loop
    already holds, and no vertical operator is applied a second time.
"""
struct ImplicitVertical{MX,MY} <: AbstractVerticalTreatment
    thomas_x::MX
    thomas_y::MY
end
Adapt.@adapt_structure ImplicitVertical

function ImplicitVertical(grid::StaggeredGrid; halo = 1)
    (; arch) = grid
    g = grid.grid
    T = eltype(g)
    return ImplicitVertical(
        _field(arch, g, NODE_ACX, T, halo),
        _field(arch, g, NODE_ACY, T, halo),
    )
end

"""
$(TYPEDSIGNATURES)

Reject a vertical treatment the grid the solver is being built on cannot serve, at
construction rather than as a silently ignored setting at solve time. Only the
depth-averaged [`PseudoTransientSolver`](@ref) constructor constrains this.
"""
_check_vertical_treatment(::ExplicitVertical) = nothing

_check_vertical_treatment(::ImplicitVertical) = throw(
    ArgumentError(
        "ImplicitVertical() is a MomentumBalance3D strategy and the depth-averaged " *
        "PseudoTransientSolver(grid) would never read it. SSA has no vertical operator at all " *
        "and DIVA integrates its own out analytically (the F₁/F₂ closure), so there is nothing " *
        "here to treat implicitly. Build the solver with " *
        "`PseudoTransientSolver(grid, BlatterPattynMomentumBalance(); " *
        "vertical_treatment = ImplicitVertical(grid))` instead.",
    ),
)

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
   (1D, 2D, 3D stencils). Read only by [`ViscosityPseudoTimeStep`](@ref), and only `ndim2`
   of the three.
 - `muB`: bulk-to-shear viscosity ratio entering the PT time step. Read only by
   [`ViscosityPseudoTimeStep`](@ref) — see [`GershgorinPseudoTimeStep`](@ref) for why the
   default `1e2` is ~101× too conservative on a depth-integrated balance.
 - `theta_v`: relaxation weight of the velocity update. Read only under
   [`FixedTuning`](@ref).
 - `gamma`: damping coefficient of the pseudo-transient rate (Sandip et al. 2024,
   Eq. 12–14; Frankel 1950). The rate is accumulated as
   `dv_new = (1 - gamma) * dv_old + r(u)`, so `gamma = 1` discards all memory and recovers
   the plain (undamped, first-order) iteration; `0 < gamma < 1` turns it into a damped
   wave, which is what buys sub-quadratic iteration-count scaling. Read only under
   [`FixedTuning`](@ref), where a good value needs the same manual scan as `theta_v`;
   [`AutotunedDynamicRelaxation`](@ref) derives it instead.
 - `abstol`: convergence tolerance on whatever `convergence` measures — a velocity change
   per iteration (m yr⁻¹) for [`VelocityIncrement`](@ref), a dimensionless residual for
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
 - `dtau_x`, `dtau_y`: pseudo-time step, per grid point, at the same node class as the
   velocity component it advances (filled by [`pseudo_dt!`](@ref)).
 - `residual_x`, `residual_y`: the raw (undamped) pseudo-transient rate `r(u)/(ρH)`,
   written unconditionally by [`dotvel!`](@ref) regardless of `gamma` — a genuinely
   separate buffer from `velocity_x_dt`/`velocity_y_dt` because damping's whole point is
   for the accumulator to *not* equal the raw rate, so the raw rate needs its own home to
   remain readable. [`pseudo_transient!`](@ref) reduces these into the `residual` it
   returns.
 - `pseudo_timestep`: an [`AbstractPseudoTimeStep`](@ref), default
   [`ViscosityPseudoTimeStep`](@ref). Selects how `dtau_x`/`dtau_y` are filled.
 - `convergence`: an [`AbstractPTConvergence`](@ref), default [`VelocityIncrement`](@ref).
   Selects what `abstol` is compared against — **and therefore what `abstol` means**.
 - `viscosity_continuation`: an [`AbstractViscosityContinuation`](@ref), default
   [`NoViscosityContinuation`](@ref).
 - `tuning`: an [`AbstractPTTuning`](@ref), default [`FixedTuning`](@ref). Selects whether
   `Δτ` and the damping are the hand-set `theta_v`/`gamma`/`pseudo_timestep` values or are
   derived from spectral estimates during the solve
   ([`AutotunedDynamicRelaxation`](@ref), which requires
   `pseudo_timestep::GershgorinPseudoTimeStep` and then ignores `theta_v`/`gamma`).
 - `friction_update`: an [`AbstractFrictionUpdate`](@ref), default
   [`ActiveFrictionUpdate`](@ref). Set to [`NoFrictionUpdate`](@ref) to bypass the basal
   friction law and hold `stress.base_x`/`base_y` fixed at whatever was written there
   before the solve (e.g. a prescribed basal stress field).
 - `vertical_treatment`: an [`AbstractVerticalTreatment`](@ref), default
   [`ExplicitVertical`](@ref). Read only on the [`MomentumBalance3D`](@ref) path; selects
   whether `∂z(µ ∂z u)` is advanced explicitly with everything else or by a per-column
   implicit line solve ([`ImplicitVertical`](@ref)).

!!! note "Two type parameters for the work arrays, not one"
    On a [`StaggeredGrid`](@ref) the x- and y-velocity work arrays are Chmy `Field`s at
    different locations (`acx`/`acy`), which are different concrete types (location is
    part of a `Field`'s type), so a single `M` shared by all four work arrays would reject
    that combination outright — hence the `MX`/`MY` split below.
"""
struct PseudoTransientSolver{
    T<:AbstractFloat,
    MX,
    MY,
    PT<:AbstractPseudoTimeStep,
    CV<:AbstractPTConvergence,
    VC<:AbstractViscosityContinuation,
    FU<:AbstractFrictionUpdate,
    TU<:AbstractPTTuning,
    DU<:AbstractDIVUpdate,
    VT<:AbstractVerticalTreatment,
} <: AbstractMomentumSolver
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
    tuning::TU
    div_update::DU
    vertical_treatment::VT
end
Adapt.@adapt_structure PseudoTransientSolver

"""
$(TYPEDSIGNATURES)

Build a [`PseudoTransientSolver`](@ref) for [`pseudo_transient!`](@ref). All work arrays
are Chmy `Field`s at `acx`/`acy` on `grid.grid2d`, matching the velocity components they
mirror (`mech.velocity.x`/`y`).

No `nz == 1` requirement here: every work array is built on `grid.grid2d`, and the unknown
the solver iterates (`velocity.depthaverage_x`/`y`) is depth-integrated for SSA and DIVA
alike. Whether a given momentum balance tolerates the grid is [`_check_momentum_grid`](@ref)'s
job, at the `pseudo_transient!` call that knows which balance it is.
"""
function PseudoTransientSolver(
    grid::StaggeredGrid;
    abstol = 1e-8,
    maxiter = 100,
    ncheck = 1,
    printout_every = typemax(Int),
    dtau_scaling = 1,
    halo = 1,
    pseudo_timestep::AbstractPseudoTimeStep = GershgorinPseudoTimeStep(),
    convergence::AbstractPTConvergence = VelocityIncrement(),
    viscosity_continuation::AbstractViscosityContinuation = NoViscosityContinuation(),
    friction_update::AbstractFrictionUpdate = ActiveFrictionUpdate(),
    tuning::AbstractPTTuning = FixedTuning(),
    div_update::AbstractDIVUpdate = NoDIVUpdate(),
    vertical_treatment::AbstractVerticalTreatment = ExplicitVertical(),
)
    _check_tuning(tuning, pseudo_timestep)
    _check_vertical_treatment(vertical_treatment)
    (; arch) = grid
    g = grid.grid2d
    T = eltype(g)
    acx() = _field(arch, g, NODE_ACX, T, halo)
    acy() = _field(arch, g, NODE_ACY, T, halo)
    return PseudoTransientSolver(
        T(abstol),
        Int(maxiter),
        Int(ncheck),
        Int(printout_every),
        T(dtau_scaling),
        acx(),
        acy(),
        acx(),
        acy(),
        acx(),
        acy(),
        acx(),
        acy(),
        pseudo_timestep,
        convergence,
        viscosity_continuation,
        friction_update,
        tuning,
        div_update,
        vertical_treatment,
    )
end

"""
$(TYPEDSIGNATURES)

Build a [`PseudoTransientSolver`](@ref) for a [`MomentumBalance3D`](@ref) solve
([`pseudo_transient!`](@ref)). All work arrays are Chmy `Field`s at `acx`/`acy` on
`grid.grid` (`ACX3`/`ACY3`), matching the column velocity components they mirror
(`mech.velocity.x`/`y`) — the only difference from the depth-averaged constructor above is
which of `grid.grid2d`/`grid.grid` the work arrays are built on. A separate method rather
than a keyword on the existing constructor (`roadmaps/blatter-pattyn.md`, Phase 1: "prefer
dispatch — the balance already decides the grid via `_check_momentum_grid`, and a keyword
lets the two disagree"), and kept as a fully independent method body so the depth-averaged
constructor above is untouched by this one's existence.

No `nz > 1` requirement here either, for the same reason as the 2D constructor: whether
`momentum` tolerates the grid is [`_check_momentum_grid`](@ref)'s job, at the
`pseudo_transient!` call that knows the `Runtime`.
"""
function PseudoTransientSolver(
    grid::StaggeredGrid,
    momentum::MomentumBalance3D;
    abstol = 1e-8,
    maxiter = 100,
    ncheck = 1,
    printout_every = typemax(Int),
    dtau_scaling = 1,
    halo = 1,
    pseudo_timestep::AbstractPseudoTimeStep = GershgorinPseudoTimeStep(),
    convergence::AbstractPTConvergence = VelocityIncrement(),
    viscosity_continuation::AbstractViscosityContinuation = NoViscosityContinuation(),
    friction_update::AbstractFrictionUpdate = ActiveFrictionUpdate(),
    tuning::AbstractPTTuning = FixedTuning(),
    div_update::AbstractDIVUpdate = NoDIVUpdate(),
    vertical_treatment::AbstractVerticalTreatment = ExplicitVertical(),
)
    _check_tuning(tuning, pseudo_timestep)
    (; arch) = grid
    g = grid.grid
    T = eltype(g)
    acx() = _field(arch, g, NODE_ACX, T, halo)
    acy() = _field(arch, g, NODE_ACY, T, halo)
    _check_vertical_scratch(vertical_treatment, acx(), acy())
    return PseudoTransientSolver(
        T(abstol),
        Int(maxiter),
        Int(ncheck),
        Int(printout_every),
        T(dtau_scaling),
        acx(),
        acy(),
        acx(),
        acy(),
        acx(),
        acy(),
        acx(),
        acy(),
        pseudo_timestep,
        convergence,
        viscosity_continuation,
        friction_update,
        tuning,
        div_update,
        vertical_treatment,
    )
end

"""
$(TYPEDSIGNATURES)

Reject an [`ImplicitVertical`](@ref) built against a *different* grid than the solver it is
handed to — which would not error on its own, it would silently write the Thomas
coefficients into the wrong shape. Checked here, where both are in hand, rather than left
to a bounds error deep inside a kernel.
"""
_check_vertical_scratch(::ExplicitVertical, acx, acy) = nothing

function _check_vertical_scratch(vt::ImplicitVertical, acx, acy)
    (size(vt.thomas_x) == size(acx) && size(vt.thomas_y) == size(acy)) || throw(
        ArgumentError(
            "ImplicitVertical's scratch fields are $(size(vt.thomas_x))/$(size(vt.thomas_y)), " *
            "but this solver's work arrays are $(size(acx))/$(size(acy)). Build the strategy " *
            "from the same grid (and `halo`) as the solver: " *
            "`ImplicitVertical(grid)`.",
        ),
    )
    eltype(vt.thomas_x) === eltype(acx) || throw(
        ArgumentError(
            "ImplicitVertical's scratch is $(eltype(vt.thomas_x)) but this solver is " *
            "$(eltype(acx)). Build the strategy from the same grid as the solver.",
        ),
    )
    return nothing
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
    DYN<:AbstractMomentumBalance,
    T<:AbstractFloat,
    VT<:AbstractVector,        # float vector type (u, u0, b)
    MT,                            # sparse matrix type (SparseMatrixCSC or CuSparseMatrix)
    PI<:AbstractVector{Int},   # perm index vector type
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

# @dev TODO: not yet implemented.
function calc_F_integral(visc_eff, H_ice, f_ice, zeta_aa, n)
    error("calc_F_integral is not yet implemented")
end

function vertically_integrated_viscosity!(N, H, μ)
    N .= H .* μ
    return nothing
end

function stagger() end

function fill_pattern!(Ai, Aj, nx, ny, i_idx, j_idx)
    k = 0
    for i = 1:nx, j = 1:ny
        im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
        nr = _ij2n_ux(i, j, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_ux(ip1, j, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_ux(i, j, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_ux(im1, j, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_ux(i, jp1, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_ux(i, jm1, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_uy(i, j, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_uy(ip1, j, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_uy(ip1, jm1, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_uy(i, jm1, nx, ny)
    end
    for i = 1:nx, j = 1:ny
        im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
        nr = _ij2n_uy(i, j, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_uy(i, jp1, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_uy(i, j, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_uy(i, jm1, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_uy(ip1, j, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_uy(im1, j, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_ux(i, jp1, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_ux(i, j, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_ux(im1, jp1, nx, ny)
        k+=1
        Ai[k]=nr
        Aj[k]=_ij2n_ux(im1, j, nx, ny)
    end
    return nothing
end

function coo_to_nzval_idx(Ai, Aj, A::SparseMatrixCSC)
    perm = Vector{Int}(undef, length(Ai))
    for k in eachindex(Ai)
        c = Aj[k]
        r = Ai[k]
        lo = A.colptr[c]
        hi = A.colptr[c+1] - 1
        perm[k] = searchsortedfirst(A.rowval, r, lo, hi, Base.Order.Forward)
    end
    return perm
end

function LinearMomentumSolver2D(
    grid::RegularGrid,
    dynamics::DYN,
    backend = CPU(),
) where {DYN<:AbstractMomentumBalance}
    T = eltype(grid.x)
    nx, ny = grid.nx, grid.ny
    dx, dy = T(grid.dx), T(grid.dy)
    dxdx_ = 1 / (dx * dx)
    dydy_ = 1 / (dy * dy)
    dxdy_ = 1 / (dx * dy)
    n_sprs = 18 * nx * ny  # 9 nonzeros/row × 2 equations (ux, uy)
    n_u = 2 * nx * ny

    i_idx = PeriodicIndexing(1, nx)
    j_idx = PeriodicIndexing(1, ny)

    # Pattern and permutation are always computed on CPU (one-time cost).
    Ai_cpu = zeros(Int, n_sprs)
    Aj_cpu = zeros(Int, n_sprs)
    fill_pattern!(Ai_cpu, Aj_cpu, nx, ny, i_idx, j_idx)
    A_cpu = sparse(Ai_cpu, Aj_cpu, ones(T, n_sprs), n_u, n_u)
    perm_cpu = coo_to_nzval_idx(Ai_cpu, Aj_cpu, A_cpu)

    # Allocate live arrays on the target backend (H→D transfer on GPU, no-op copy on CPU).
    u = KernelAbstractions.zeros(backend, T, n_u)
    u0 = KernelAbstractions.zeros(backend, T, n_u)
    b = KernelAbstractions.zeros(backend, T, n_u)
    perm = KernelAbstractions.zeros(backend, Int, n_sprs)
    perm .= perm_cpu

    # `A_cpu` stays a `SparseMatrixCSC` for the CPU default; for GPU, adapt it before
    # passing to the inner constructor, e.g. `CUDA.CUSPARSE.CuSparseMatrixCSC(A_cpu)`.
    return LinearMomentumSolver2D(
        dynamics,
        nx,
        ny,
        dxdx_,
        dydy_,
        dxdy_,
        u,
        u0,
        b,
        A_cpu,
        perm,
        i_idx,
        j_idx,
        Ref{Any}(nothing),
    )
end

@kernel function _assemble_ux!(
    nzval,
    perm,
    u0,
    b,
    nx,
    ny,
    N,
    N_ab,
    ux,
    taud_acx,
    β_acx,
    β_acy,
    dxdx_,
    dydy_,
    dxdy_,
    i_idx,
    j_idx,
    dynamics,
)
    i, j = @index(Global, NTuple)
    im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
    nr = _ij2n_ux(i, j, nx, ny)
    @inbounds u0[nr] = ux[i, j]
    @inbounds b[nr] = taud_acx[i, j]
    v = loop1_coeffs(
        im1,
        i,
        ip1,
        jm1,
        j,
        jp1,
        dxdx_,
        dydy_,
        dxdy_,
        N,
        N_ab,
        β_acx,
        β_acy,
        dynamics,
    )
    k0 = 9 * ((i - 1) * ny + (j - 1))
    @inbounds for s = 1:9
        nzval[perm[k0+s]] = v[s]
    end
end

@kernel function _assemble_uy!(
    nzval,
    perm,
    u0,
    b,
    nx,
    ny,
    N,
    N_ab,
    uy,
    taud_acy,
    β_acx,
    β_acy,
    dxdx_,
    dydy_,
    dxdy_,
    i_idx,
    j_idx,
    dynamics,
)
    i, j = @index(Global, NTuple)
    im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
    nr = _ij2n_uy(i, j, nx, ny)
    @inbounds u0[nr] = uy[i, j]
    @inbounds b[nr] = taud_acy[i, j]
    v = loop2_coeffs(
        im1,
        i,
        ip1,
        jm1,
        j,
        jp1,
        dxdx_,
        dydy_,
        dxdy_,
        N,
        N_ab,
        β_acx,
        β_acy,
        dynamics,
    )
    k0 = 9 * nx * ny + 9 * ((i - 1) * ny + (j - 1))
    @inbounds for s = 1:9
        nzval[perm[k0+s]] = v[s]
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
    N,
    N_ab,
    ux,
    uy,
    taud_acx,
    taud_acy,
    β_acx,
    β_acy,
)
    (; A, perm, u0, b, i_idx, j_idx, nx, ny, dxdx_, dydy_, dxdy_, dynamics) = lsd
    backend = get_backend(lsd.u)
    kern_ux = _assemble_ux!(backend)
    kern_uy = _assemble_uy!(backend)

    nzval = nonzeros(A)
    kern_ux(
        nzval,
        perm,
        u0,
        b,
        nx,
        ny,
        N,
        N_ab,
        ux,
        taud_acx,
        β_acx,
        β_acy,
        dxdx_,
        dydy_,
        dxdy_,
        i_idx,
        j_idx,
        dynamics;
        ndrange = (nx, ny),
    )
    kern_uy(
        nzval,
        perm,
        u0,
        b,
        nx,
        ny,
        N,
        N_ab,
        uy,
        taud_acy,
        β_acx,
        β_acy,
        dxdx_,
        dydy_,
        dxdy_,
        i_idx,
        j_idx,
        dynamics;
        ndrange = (nx, ny),
    )
    return nothing
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
    @inbounds for i = 1:nx, j = 1:ny
        ux[i, j] = u[_ij2n_ux(i, j, nx, ny)]
        uy[i, j] = u[_ij2n_uy(i, j, nx, ny)]
    end
    return nothing
end
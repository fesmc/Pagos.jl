###############################################################
# Pseudo-transient momentum solver (Sandip et al., 2024)
###############################################################

"""
$(TYPEDSIGNATURES)

Relax the velocity field: `v = v_old + theta_v * dotvel * dtau`.
"""
function pseudo_vel!(v, v_old, dotvel, dtau, theta_v)
    @. v = v_old + theta_v * dotvel * dtau
    return nothing
end

###############################################################
# C-grid staggered pseudo-transient solver
###############################################################
#
# Two facts govern everything below and neither lives inside a kernel. `Chmy.Field` is an
# `AbstractArray` with no `BroadcastStyle`, so `copyto!`, the `v = v_old + θ·dv·dτ`
# relaxation and the max-norm reductions would fall back to scalar `getindex` — right on
# CPU, catastrophic on GPU — and are therefore all routed through `asarray`. And the
# membrane stress reads velocity gradients one ring beyond the interior at `ab`, so the
# velocity halo must be refreshed every iteration, not once at the end; `Neumann(0)` is the
# placeholder until per-equation boundary conditions land (`roadmaps/chmy.md` Phase 4).

"""
$(TYPEDSIGNATURES)

Local (per-grid-point) pseudo-transient time step (Sandip et al. 2024, Eq. 7): writes
`dtau_x`, `dtau_y` from the depth-averaged viscosity `mu` (`aa`), staggered onto the
velocity's own node class (`acx`/`acy`) by `lerp` — the same arithmetic-mean choice
[`drivingstress!`](@ref) and [`dotvel!`](@ref) make for the ice thickness `H`, and folding
the safety factor `dtau_scaling` in up front.

Per grid point rather than a `maximum(mu)` reduction, so a single stiff cell no longer
throttles `Δτ` everywhere else and no global reduction is needed. Only `mu`'s interior is
read (`lerp` needs its halo only at the two domain-boundary `ab`-style corners the membrane
stress itself needs — see [`pseudo_transient!`](@ref)'s halo warning), so this call has the
same halo requirement as the rest of the solve.
"""
function pseudo_dt!(dtau_x, dtau_y, ρ, dx, dy, mu, muB, ndim, dtau_scaling, rt::Runtime)
    scaling = convert(eltype(dtau_x), dtau_scaling * ρ * dx * dy / (4 * (1 + muB) * ndim))
    rt.launch2d(rt.arch, rt.grid2d,
              _pseudo_dt_staggered! => (dtau_x, dtau_y, mu, scaling, rt.grid2d))
    return nothing
end

@kernel inbounds = true function _pseudo_dt_staggered!(dtau_x, dtau_y, mu, scaling, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    dtau_x[I...] = scaling / lerp(mu, NODE_ACX, grid, I...)
    dtau_y[I...] = scaling / lerp(mu, NODE_ACY, grid, I...)
end

###############################################################
# Gershgorin pseudo-time step
###############################################################
#
# `ViscosityPseudoTimeStep` (Sandip Eq. 7) bounds the *membrane* part of the operator and
# nothing else. That is fine on the uniform-slab tests, where β is small and there is no ice
# margin, and wrong in both directions on real geometry:
#
#  - it omits the basal-drag term β/(ρH) from the bound, which under grounded Antarctic ice
#    (β_eff up to ~5e13 Pa s m⁻¹) is the *dominant* eigenvalue — a solve that recomputes
#    friction from the velocity (`ActiveFrictionUpdate`) diverges within ~20 iterations;
#  - it `lerp`s the viscosity straight across the ice margin, so an ice-free cell's
#    placeholder viscosity throttles Δτ on exactly the faces that carry the calving-front
#    forcing.
#
# The bound below is Gershgorin's: λ_max ≤ max row sum of |coefficients| of the mass-scaled
# residual operator, and explicit stability is Δτ ≤ 2/λ_max. The coefficients are read off
# the discretisation the solver actually runs — `_membrane_stress_staggered!`'s `ηH` at `aa`
# and `hlerp(η)·lerp(H)` at `ab`, `_basalstress_staggered!`'s `lerp(β)`, and
# `_dotvel_staggered!`'s `lerp(H)` — and every one of them is evaluated through the same
# mask those kernels use, so an ice-free neighbour contributes exactly the zero it
# contributes to the residual itself.
#
# Row sum for the u-equation at `acx(i, j)`, writing P = ηH at `aa` and Q = hlerp(η)·lerp(H)
# at `ab` (derivation: expand ∂x(N_xx) + ∂y(N_xy) - βu and sum |coefficients| over both the
# u and the v unknowns, since the two equations are coupled):
#
#   ∂x(N_xx) : u-terms 8(P₋+P₊)/dx²   v-terms 4(P₋+P₊)/(dx·dy)
#   ∂y(N_xy) : u-terms 2(Q₋+Q₊)/dy²   v-terms 2(Q₋+Q₊)/(dx·dy)
#   basal    : β_face
#
# all divided by ρ·H_face. Sanity check: uniform η, H, no drag, dx = dy gives
# Λ = 32η/(ρdx²) hence Δτ = ρdx²/(16η) — which is Sandip's Eq. 7 at muB = 0, ndim2 = 4.1
# (ρdx²/16.4η). The two agree where they should; they differ only by the terms Eq. 7 omits.

@inline _etaH_aa(η, H, mask, i, j, k) =
    node_active(mask, NODE_AA, i, j) ? η[i, j, k] * H[i, j, k] : zero(eltype(η))

# `node_fully_active`, matching `_membrane_stress_staggered!`: `hlerp` of a zero-viscosity
# neighbour is `NaN`, not zero, so the strict rule is what makes this evaluable at all.
@inline _etaH_ab(η, H, mask, grid, i, j, k) =
    node_fully_active(mask, NODE_AB, i, j) ?
    hlerp(η, NODE_AB, grid, i, j, k) * lerp(H, NODE_AB, grid, i, j, k) : zero(eltype(η))

# Drag enters the bound only if the drag term actually depends on the velocity being
# iterated. Under `NoFrictionUpdate` the basal stress is a prescribed constant field, i.e. a
# forcing like the driving stress, contributing nothing to the operator's spectrum.
@inline _drag_in_spectrum(::ActiveFrictionUpdate) = true
@inline _drag_in_spectrum(::NoFrictionUpdate) = false

@kernel inbounds = true function _pseudo_dt_gershgorin!(dtau_x, dtau_y, η, H, β, ρ, scale,
                                                        drag, mask, dx, dy, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, k = I
    Z = zero(eltype(dtau_x))

    if node_active(mask, NODE_ACX, i, j)
        Hx = lerp(H, NODE_ACX, grid, I...)
        Ps = _etaH_aa(η, H, mask, i - 1, j, k) + _etaH_aa(η, H, mask, i, j, k)
        Qs = _etaH_ab(η, H, mask, grid, i, j, k) + _etaH_ab(η, H, mask, grid, i, j + 1, k)
        drag_x = drag ? lerp(β, NODE_ACX, grid, I...) : Z
        Λ = (8Ps / dx^2 + 4Ps / (dx * dy) + 2Qs / dy^2 + 2Qs / (dx * dy) + drag_x) / (ρ * Hx)
        dtau_x[I...] = (Hx > Z && Λ > Z) ? scale / Λ : Z
    else
        dtau_x[I...] = Z
    end

    if node_active(mask, NODE_ACY, i, j)
        Hy = lerp(H, NODE_ACY, grid, I...)
        Ps = _etaH_aa(η, H, mask, i, j - 1, k) + _etaH_aa(η, H, mask, i, j, k)
        Qs = _etaH_ab(η, H, mask, grid, i, j, k) + _etaH_ab(η, H, mask, grid, i + 1, j, k)
        drag_y = drag ? lerp(β, NODE_ACY, grid, I...) : Z
        Λ = (8Ps / dy^2 + 4Ps / (dx * dy) + 2Qs / dx^2 + 2Qs / (dx * dy) + drag_y) / (ρ * Hy)
        dtau_y[I...] = (Hy > Z && Λ > Z) ? scale / Λ : Z
    else
        dtau_y[I...] = Z
    end
end

"""
$(TYPEDSIGNATURES)

Fill `solver.dtau_x`/`dtau_y` for one solve, dispatching on `solver.pseudo_timestep`:
[`ViscosityPseudoTimeStep`](@ref) (the default) delegates to the Sandip Eq. 7 method above,
[`GershgorinPseudoTimeStep`](@ref) bounds the spectral radius of the residual operator the
solver actually iterates — including the basal-drag term and the ice mask (see the source
note above for the row sums, and [`GershgorinPseudoTimeStep`](@ref) for why they matter).

Called once by [`pseudo_transient!`](@ref) before the PT loop; `mask` and
`solver.friction_update` must be the same ones the loop itself will use, or the bound
describes a different operator than the one being iterated.
"""
function pseudo_dt!(solver::PseudoTransientSolver, mech::MechanicState, c::Constants,
                    rt::Runtime, mask::AbstractIceMask = NoMask())
    dx = Δx(rt.grid2d, Center(), 1, 1, 1)
    dy = Δy(rt.grid2d, Center(), 1, 1, 1)
    return pseudo_dt!(solver, solver.pseudo_timestep, mech, c, rt, mask, dx, dy)
end

function pseudo_dt!(solver::PseudoTransientSolver, pt::ViscosityPseudoTimeStep,
                    mech::MechanicState, c::Constants, rt::Runtime,
                    ::AbstractIceMask, dx, dy)
    pseudo_dt!(solver.dtau_x, solver.dtau_y, c.density_ice, dx, dy,
               mech.material.viscosity_depthaveraged, pt.muB, pt.ndim,
               solver.dtau_scaling, rt)
    return nothing
end

function pseudo_dt!(solver::PseudoTransientSolver, pt::GershgorinPseudoTimeStep,
                    mech::MechanicState, c::Constants, rt::Runtime,
                    mask::AbstractIceMask, dx, dy)
    gershgorin_dt!(solver, mech, c, rt, mask, 2 * pt.cfl * solver.dtau_scaling)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Fill `solver.dtau_x`/`dtau_y` with `scale / Λ`, `Λ` being the Gershgorin absolute row sum of
the mass-scaled residual operator (see the source note above). The prefactor is a parameter
because the two schemes reading this bound want different ones: the first-order iteration
takes `scale = 2·cfl` (`Δτ ≤ 2/λ_max`), [`AutotunedDynamicRelaxation`](@ref) takes
`scale = Δτ²` (`Δτ` enters the damped update twice, on the residual and on the accumulator).

`Λ = scale / dtau` is recoverable from the output, which is how the autotuner gets its
`M`-inner product without a third field.
"""
function gershgorin_dt!(solver::PseudoTransientSolver, mech::MechanicState, c::Constants,
                        rt::Runtime, mask::AbstractIceMask, scale)
    T = eltype(solver.dtau_x)
    dx = Δx(rt.grid2d, Center(), 1, 1, 1)
    dy = Δy(rt.grid2d, Center(), 1, 1, 1)
    rt.launch2d(rt.arch, rt.grid2d,
                _pseudo_dt_gershgorin! =>
                    (solver.dtau_x, solver.dtau_y, mech.material.viscosity_depthaveraged,
                     mech.topography.thickness, mech.friction.beta_eff,
                     convert(T, c.density_ice), convert(T, scale),
                     _drag_in_spectrum(solver.friction_update), mask,
                     convert(T, dx), convert(T, dy), rt.grid2d))
    return nothing
end

"""
$(TYPEDSIGNATURES)

The SSA/DIVA pseudo-transient velocity rate. The membrane-stress divergence is formed by Chmy's `∂x`/`∂y` acting directly on the
already-staggered `sxx` (`aa`), `sxy` (`ab`), `syy` (`aa`) — no interpolation, by the same
node-algebra argument as [`drivingstress!`](@ref) and [`velocitygradients!`](@ref): `∂x`
of `sxx` and `∂y` of `sxy` both land on `acx`; `∂x` of `sxy` and `∂y` of `syy` both land on
`acy`. Only the ice thickness `H` needs staggering onto the face, by `lerp`.

Depth-integrated throughout, so it runs on `rt.grid2d`.

The rate is accumulated with damping (Sandip et al. 2024, Eq. 12–14):
`dvx = (1 - gamma) * dvx_old + r_x(u)`, where `dvx_old` is whatever `dvx` already holds on
entry. `gamma = 1` (the default) discards `dvx_old` entirely and reproduces the plain,
undamped rate every call — this method's behaviour before `gamma` was added.
`0 < gamma < 1` gives the rate memory, which is what
lets [`pseudo_transient!`](@ref) converge in fewer iterations (the whole point of Sandip's
"acceleration owing to damping", their §2.5).

`resid_x`, `resid_y` are mandatory keyword outputs (no default: aliasing `dvx`/`dvy` would
corrupt the damped combination above, so a genuinely separate buffer is required) that
always receive the raw, undamped rate `r(u)/(ρH)` — independent of `gamma`, and thus a
meaningful "how far from momentum balance" diagnostic even when the accumulator itself is
damped. [`pseudo_transient!`](@ref) reduces these into the `residual` it returns.
"""
function dotvel!(dvx, dvy, sxx, sxy, syy, base_x, base_y, driving_x, driving_y,
                 H, density_ice, rt::Runtime,
                 momentum::MomentumBalance2D,
                 mask::AbstractIceMask = NoMask(); gamma = 1, resid_x, resid_y)
    ρ = convert(eltype(dvx), density_ice)
    γ = convert(eltype(dvx), gamma)
    rt.launch2d(rt.arch, rt.grid2d,
              _dotvel_staggered! => (dvx, dvy, resid_x, resid_y, sxx, sxy, syy, base_x,
                                     base_y, driving_x, driving_y, H, ρ, γ, mask, rt.grid2d))
    return nothing
end

@kernel inbounds = true function _dotvel_staggered!(dvx, dvy, resid_x, resid_y, sxx, sxy,
                                                     syy, base_x, base_y, driving_x,
                                                     driving_y, H, ρ, γ, mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    Z = zero(eltype(dvx))
    one_m_γ = one(γ) - γ

    if node_active(mask, NODE_ACX, i, j)
        Hx = lerp(H, NODE_ACX, grid, I...)
        shear_x = ∂x(sxx, grid, I...) + ∂y(sxy, grid, I...)
        raw_x = Hx > zero(Hx) ?
                (shear_x - base_x[I...] - driving_x[I...]) / (ρ * Hx) : Z
        resid_x[I...] = raw_x
        dvx[I...] = one_m_γ * dvx[I...] + raw_x
    else
        resid_x[I...] = Z
        dvx[I...] = Z
    end

    if node_active(mask, NODE_ACY, i, j)
        Hy = lerp(H, NODE_ACY, grid, I...)
        shear_y = ∂x(sxy, grid, I...) + ∂y(syy, grid, I...)
        raw_y = Hy > zero(Hy) ?
                (shear_y - base_y[I...] - driving_y[I...]) / (ρ * Hy) : Z
        resid_y[I...] = raw_y
        dvy[I...] = one_m_γ * dvy[I...] + raw_y
    else
        resid_y[I...] = Z
        dvy[I...] = Z
    end
end

###############################################################
# Convergence criteria
###############################################################

# `_pt_error` is what the PT loop compares against `abstol`. Both methods reduce fields the
# loop already maintains, so neither costs an extra kernel — only the reduction the
# `ncheck` cadence exists to amortize.

@inline _abs_diff(a, b) = abs(a - b)

_pt_error(::VelocityIncrement, solver, ux, uy, ref) =
    max(mapreduce(_abs_diff, max, asarray(ux), asarray(solver.velocity_x_old)),
        mapreduce(_abs_diff, max, asarray(uy), asarray(solver.velocity_y_old)))

_pt_error(::ScaledResidual, solver, ux, uy, ref) =
    max(maximum(abs, asarray(solver.residual_x)),
        maximum(abs, asarray(solver.residual_y))) / ref

@kernel inbounds = true function _driving_rate!(rate_x, rate_y, driving_x, driving_y, H, ρ,
                                                 mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    Z = zero(eltype(rate_x))
    Hx = lerp(H, NODE_ACX, grid, I...)
    Hy = lerp(H, NODE_ACY, grid, I...)
    rate_x[I...] = (node_active(mask, NODE_ACX, i, j) && Hx > Z) ?
                   driving_x[I...] / (ρ * Hx) : Z
    rate_y[I...] = (node_active(mask, NODE_ACY, i, j) && Hy > Z) ?
                   driving_y[I...] / (ρ * Hy) : Z
end

# The scale `_pt_error(::ScaledResidual, ...)` divides by: the velocity rate the driving
# stress alone would produce, i.e. `dotvel!`'s residual with the membrane and basal terms
# dropped. Computed once per solve, before the loop.
#
# `solver.residual_x`/`residual_y` are borrowed as scratch: `dotvel!` overwrites them
# unconditionally on the first iteration, so nothing that is read afterwards is lost, and
# the alternative — two more `acx`/`acy` fields on every solver ever constructed — would be
# paid for by every solve to serve one convergence criterion.
_convergence_scale(::VelocityIncrement, mech, c, solver, rt, mask) =
    one(eltype(asarray(mech.velocity.depthaverage_x)))

function _convergence_scale(::ScaledResidual, mech::MechanicState, c::Constants,
                            solver::PseudoTransientSolver, rt::Runtime,
                            mask::AbstractIceMask)
    T = eltype(solver.residual_x)
    rt.launch2d(rt.arch, rt.grid2d,
                _driving_rate! => (solver.residual_x, solver.residual_y,
                                   mech.stress.driving_x, mech.stress.driving_y,
                                   mech.topography.thickness, convert(T, c.density_ice),
                                   mask, rt.grid2d))
    scale = max(maximum(abs, asarray(solver.residual_x)),
                maximum(abs, asarray(solver.residual_y)))
    # A domain with no driving stress at all (a flat, ice-free or perfectly level state) is
    # already in balance; falling back to 1 makes `err` the raw residual rather than `Inf`.
    return scale > zero(scale) ? scale : one(scale)
end

###############################################################
# Autotuned dynamic relaxation (Duretz et al. 2026)
###############################################################

"""
$(TYPEDSIGNATURES)

The iteration parameters in force at one point of the PT loop, as an immutable
`NamedTuple`: `gamma`, `theta_v`, the `Δτ` prefactor `scale` (`dtau = scale/Λ`, see
[`gershgorin_dt!`](@ref)), the last `lambda_min`, and the two partial sums
[`_arm_tuning`](@ref)/[`_tune!`](@ref) carry between iterations.

The field set is the same for every [`AbstractPTTuning`](@ref), so [`pseudo_transient!`](@ref)
stays type stable whichever is selected. `scale`/`lambda_min` are `NaN` under
[`FixedTuning`](@ref), which reads neither.
"""
_tuning_state(::Type{T}, gamma, theta_v, scale) where {T} =
    (; gamma = T(gamma), theta_v = T(theta_v), scale = T(scale), lambda_min = T(NaN),
       rayleigh_ur = zero(T), rayleigh_uu = zero(T), armed = false)

"""
$(TYPEDSIGNATURES)

Fill `solver.dtau_x`/`dtau_y` for the first iteration and return the initial
[`_tuning_state`](@ref).

[`FixedTuning`](@ref) defers to [`pseudo_dt!`](@ref) and freezes its own `gamma`/`theta_v`.
[`AutotunedDynamicRelaxation`](@ref) starts in its warm-up: undamped (`γ = 1`) at the
damped-Jacobi step `Δτ = 1/λ_max`, i.e. `scale = 1`.
"""
function _tuning_init!(tu::FixedTuning, solver::PseudoTransientSolver, mech::MechanicState,
                       c::Constants, rt::Runtime, mask::AbstractIceMask)
    pseudo_dt!(solver, mech, c, rt, mask)
    return _tuning_state(eltype(solver.dtau_x), tu.gamma, tu.theta_v, NaN)
end

function _tuning_init!(::AutotunedDynamicRelaxation, solver::PseudoTransientSolver,
                       mech::MechanicState, c::Constants, rt::Runtime,
                       mask::AbstractIceMask)
    T = eltype(solver.dtau_x)
    scale = T(solver.dtau_scaling)
    gershgorin_dt!(solver, mech, c, rt, mask, scale)
    return _tuning_state(T, 1, 1, scale)
end

@inline _du_dot_r(u, u_old, r) = (u - u_old) * r

# Guarded: `dtau` is exactly zero off-mask, where `Δu` is zero too, and an unguarded 0/0
# would `NaN` the whole reduction.
@inline _du2_over_dtau(u, u_old, dtau) =
    dtau > zero(dtau) ? (u - u_old)^2 / dtau : zero(dtau)

"""
$(TYPEDSIGNATURES)

`Δuᵀ r̃` over both velocity components, with `Δu = ux - solver.velocity_x_old`: the
numerator of the Rayleigh quotient, sampled at the two iterates that
[`_arm_tuning`](@ref) and [`_tune!`](@ref) each see one of.
"""
_sum_du_dot_r(solver, ux, uy) =
    mapreduce(_du_dot_r, +, asarray(ux), asarray(solver.velocity_x_old),
              asarray(solver.residual_x)) +
    mapreduce(_du_dot_r, +, asarray(uy), asarray(solver.velocity_y_old),
              asarray(solver.residual_y))

"""
$(TYPEDSIGNATURES)

`Σ Δu²/Δτ`, the Rayleigh quotient's `M`-inner product `Δuᵀ M Δu` up to the factor `scale`:
`M = diag(Λ)` is never stored, but [`gershgorin_dt!`](@ref) writes `dtau = scale/Λ`, so
`Λ = scale/dtau`.
"""
_sum_du2_over_dtau(solver, ux, uy) =
    mapreduce(_du2_over_dtau, +, asarray(ux), asarray(solver.velocity_x_old),
              asarray(solver.dtau_x)) +
    mapreduce(_du2_over_dtau, +, asarray(uy), asarray(solver.velocity_y_old),
              asarray(solver.dtau_y))

"""
$(TYPEDSIGNATURES)

Sample `Δuᵀ r̃(u^{k-1})` and `Δuᵀ M Δu`, and arm [`_tune!`](@ref) to take the missing
`Δuᵀ r̃(u^k)` next iteration. Splitting the Rayleigh quotient this way is what lets it be
evaluated from iterates the loop already holds, carrying two scalars instead of a second
residual buffer on every solver ever constructed.

Call after `pseudo_vel!`, where `ux`/`uy` hold `u^k`, `solver.velocity_*_old` hold
`u^{k-1}` and `solver.residual_*` hold `r̃(u^{k-1})`. A no-op under [`FixedTuning`](@ref),
and except every `cadence` iterations.
"""
_arm_tuning(::FixedTuning, state, solver, ux, uy, iter) = state

function _arm_tuning(tu::AutotunedDynamicRelaxation, state, solver::PseudoTransientSolver,
                     ux, uy, iter::Int)
    iter % tu.cadence == 0 || return state
    return merge(state, (; rayleigh_ur = _sum_du_dot_r(solver, ux, uy),
                           rayleigh_uu = _sum_du2_over_dtau(solver, ux, uy),
                           armed = true))
end

"""
$(TYPEDSIGNATURES)

Close the Rayleigh quotient armed by [`_arm_tuning`](@ref), derive `Δτ` and `γ` from it
(closed form in [`AutotunedDynamicRelaxation`](@ref)) and refill
`solver.dtau_x`/`dtau_y` — the refill doubling as the `λ_max` re-estimation a solve with an
evolving viscosity needs, since it rebuilds the Gershgorin bound from `η` as it now stands.

Call after `pseudo_rate!` and *before* the `u → u_old` copy, the one point where `Δu` is
still the increment [`_arm_tuning`](@ref) measured against and `solver.residual_*` already
hold `r̃(u^k)`. A no-op under [`FixedTuning`](@ref), and on a degenerate quotient (`Δu = 0`
at a stationary iterate), where the previous parameters stand.
"""
_tune!(::FixedTuning, state, solver, mech, c, rt, mask, ux, uy) = state

function _tune!(tu::AutotunedDynamicRelaxation, state, solver::PseudoTransientSolver,
                mech::MechanicState, c::Constants, rt::Runtime, mask::AbstractIceMask,
                ux, uy)
    state.armed || return state
    T = typeof(state.gamma)
    numerator = abs(_sum_du_dot_r(solver, ux, uy) - state.rayleigh_ur)
    denominator = state.scale * state.rayleigh_uu
    (numerator > 0 && denominator > 0) || return merge(state, (; armed = false))

    # A Rayleigh quotient of `Â` cannot exceed `λ_max ≤ 1`; above that it is measuring the
    # roundoff a converged `Δu` degenerates into, not the spectrum.
    λ_min = min(T(numerator / denominator), one(T))
    # `convert`ed, not promoted: `tuning` and `pseudo_timestep` are independent keyword
    # arguments, and a Float64 `cfl` on a Float32 solver must not retype `state`.
    d = 2 * convert(T, tu.c_damp) * sqrt(λ_min)
    cc = convert(T, solver.pseudo_timestep.cfl)^2
    dtau = -cc * d + sqrt(cc^2 * d^2 + 4 * cc)
    scale = dtau^2 * convert(T, solver.dtau_scaling)
    gershgorin_dt!(solver, mech, c, rt, mask, scale)
    return merge(state, (; gamma = d * dtau, scale, lambda_min = λ_min, armed = false))
end

###############################################################
# Glen viscosity continuation
###############################################################
#
# The dispatch types live in `src/mechanics/solvers.jl` next to `PseudoTransientSolver`;
# `update_viscosity!` lives here because it is PT-loop behaviour, called from
# `pseudo_rate!` below.

"""
$(TYPEDSIGNATURES)

No-op: `material.viscosity_depthaveraged` is untouched, since
[`NoViscosityContinuation`](@ref) means "treat viscosity as a fixed input".
"""
update_viscosity!(mech::MechanicState, ::NoViscosityContinuation, rt::Runtime,
                  mask::AbstractIceMask = NoMask()) = nothing

"""
$(TYPEDSIGNATURES)

[`GlenViscosityContinuation`](@ref): overwrite `material.viscosity_depthaveraged` in place
from the current velocity iterate — the effective strain rate
([`effective_strainrate_ssa!`](@ref), which requires `velocity`'s gradients to already be
current, i.e. [`velocitygradients!`](@ref) must run first) and `material.rate_factor_depthaveraged`,
relaxed in log-space toward the field's own previous value (see
[`GlenViscosityContinuation`](@ref)'s docstring for the formula).
"""
function update_viscosity!(mech::MechanicState, vc::GlenViscosityContinuation, rt::Runtime,
                          mask::AbstractIceMask = NoMask())
    (; velocity, strainrate, material) = mech
    effective_strainrate_ssa!(strainrate, velocity, rt, mask)
    rt.launch2d(rt.arch, rt.grid2d,
              _glen_viscosity_continuation! =>
                  (material.viscosity_depthaveraged, material.rate_factor_depthaveraged,
                   strainrate.effective_depthaveraged, vc.n_glen, vc.strainrate_reg,
                   vc.theta_mu, mask))
    return nothing
end

# No `grid` argument: every field here is already at `aa`, so unlike this file's other
# kernels there is no `lerp`/`hlerp` staggering to do.
@kernel inbounds = true function _glen_viscosity_continuation!(μ, A, eff, n, ε̇0, θ, mask, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    if node_active(mask, NODE_AA, i, j)
        ε̇ = sqrt(eff[I...]^2 + ε̇0^2)
        μ_raw = inv(2 * A[I...]^(1 / n)) * ε̇^((1 - n) / n)
        μ[I...] = exp(θ * log(μ_raw) + (1 - θ) * log(μ[I...]))
    end
end

"""
$(TYPEDSIGNATURES)

[`DIVAViscosityContinuation`](@ref): the same Glen law as the depth-averaged method above,
evaluated **per layer** on the column grid. Writes `material.viscosity` (`µ(z)`) from
DIVA's effective strain rate ([`effective_strainrate_diva!`](@ref), Eq. 13) and the column
rate factor `material.rate_factor`, then derives `material.viscosity_depthaveraged` from it
by [`depthaverage!`](@ref).

`µ̄` is *derived*, never computed independently: the membrane stress and the vertical shear
must describe the same ice, so `µ̄` has to be the average of exactly the `µ(z)` the shear was
built from.

Reuses `_glen_viscosity_continuation!` unchanged — the kernel indexes every field at its own
`I`, so the identical code serves the `aa`/`grid2d` and `aa`/`grid` launches.
"""
function update_viscosity!(mech::MechanicState, vc::DIVAViscosityContinuation, rt::Runtime,
                           mask::AbstractIceMask = NoMask())
    (; material, strainrate) = mech
    effective_strainrate_diva!(mech, rt, mask)
    rt.launch(rt.arch, rt.grid,
              _glen_viscosity_continuation! =>
                  (material.viscosity, material.rate_factor, strainrate.effective,
                   vc.n_glen, vc.strainrate_reg, vc.theta_mu, mask))
    depthaverage!(material.viscosity_depthaveraged, material.viscosity, rt, mask)
    return nothing
end

###############################################################
# DIVA depth-integrated-viscosity chain
###############################################################
#
# The full chain, in dependency order:
#
#   µ(z)  ←  ε̇_e (Eq. 13, u_z from Eq. 21 using the *previous* τ_b and µ)
#   µ̄     ←  depthaverage(µ(z))                                   [both in update_viscosity!]
#   F₁,F₂ ←  ∫(1/µ)((s-z)/H)^m dz   (Eq. 15)
#   β_eff ←  1/(1/β + F₂)           (Eqs. 19–20)
#
# `τ_b` is deliberately *not* recomputed here: it is read as whatever the last
# `update_basalstress!` left, which is exactly the paper's "obtained from the previous
# iteration". On a cold start it is zero, so `u_z = 0` and the first `ε̇_e` is the SSA
# invariant — a sane starting point rather than a singular one.

"""
$(TYPEDSIGNATURES)

Evaluate DIVA's depth-integrated-viscosity chain once, in dependency order: `µ(z)` and `µ̄`
via [`update_viscosity!`](@ref), then `F₁`/`F₂` via [`viscosity_integrals!`](@ref), then
`β_eff` via [`beta_eff_diva!`](@ref).

This is the function a caller must invoke before [`pseudo_transient!`](@ref) under the
default [`NoDIVUpdate`](@ref) — the solver will not do it for you (`roadmaps/chmy.md`,
Phase 3, decision 7). Under [`PeriodicDIVUpdate`](@ref) the loop also calls it every
`n_update` iterations.

!!! warning "Skipping this leaves `β_eff = 0`, i.e. frictionless sliding"
    `friction.beta_eff` is zero at allocation, and a zero friction coefficient is a
    perfectly well-formed (if physically absurd) input that raises no error — the solve will
    simply run away. There is deliberately no runtime guard; this docstring is the contract,
    exactly as with the halo-filling requirement on [`pseudo_transient!`](@ref).

Requires `mech.material.viscosity` to hold a usable previous iterate (the chain's `ε̇_e`
divides by it) and the depth-averaged velocity gradients to be current.
"""
function diva_update!(mech::MechanicState, solver::PseudoTransientSolver, rt::Runtime,
                      mask::AbstractIceMask = NoMask())
    (; material) = mech
    update_viscosity!(mech, solver.viscosity_continuation, rt, mask)
    viscosity_integrals!(material.viscosity_integral_1, material.viscosity_integral_2,
                         mech, rt, mask)
    beta_eff_diva!(mech, rt, mask)
    return nothing
end

# "How often", per decision 6 — kept strictly separate from "how" (the continuation type).
# SSA has no DIV chain, so its viscosity continuation keeps its every-iteration cadence
# inside `pseudo_rate!` untouched; only the DIVA path consults `solver.div_update`.
_div_refresh_due(::NoDIVUpdate, iter) = false
_div_refresh_due(d::PeriodicDIVUpdate, iter) = iter % d.n_update == 0

# Who drives the viscosity from inside `pseudo_rate!`, i.e. every PT iteration.
#
# SSA: its continuation (`GlenViscosityContinuation`) writes `viscosity_depthaveraged`
# directly from Eq. 12, and doing so every iteration *is* the continuation — that is what
# lets the nonlinear solve relax rather than jump. Unchanged.
#
# DIVA: nothing here. Its continuation writes `viscosity` (`µ(z)`) and only `diva_update!`
# may drive it, at the cadence `solver.div_update` names. Calling it here as well would
# refresh `µ(z)` every iteration regardless of that strategy — the exact contradiction the
# "how" / "how often" split exists to prevent. The two continuations write *different*
# fields (`AA2` vs `AA3`), which is why SSA's cadence needs no gate at all.
_iterate_viscosity!(mech::MechanicState, ::SSAMomentumBalance,
                    solver::PseudoTransientSolver, rt::Runtime, mask::AbstractIceMask) =
    update_viscosity!(mech, solver.viscosity_continuation, rt, mask)

_iterate_viscosity!(::MechanicState, ::DIVAMomentumBalance, ::PseudoTransientSolver,
                    ::Runtime, ::AbstractIceMask) = nothing

###############################################################
# Basal friction update
###############################################################
#
# Split between `src/mechanics/solvers.jl` and here for the same reason as the
# viscosity-continuation trio above.

"""
$(TYPEDSIGNATURES)

[`ActiveFrictionUpdate`](@ref): the ordinary basal friction law. Overwrites
`velocity.base_x`/`base_y` from the current (SSA-limit) depth-averaged velocity, then
`stress.base_x`/`base_y` from `friction.beta_eff * velocity.base_{x,y}` via
[`basalstress!`](@ref) — the behaviour every solver had before
[`AbstractFrictionUpdate`](@ref) existed.
"""
function update_basalstress!(mech::MechanicState, ::ActiveFrictionUpdate, rt::Runtime,
                             mask::AbstractIceMask = NoMask())
    (; velocity, stress, friction) = mech
    # The SSA limit `u_b = ū`, which is what both balances use today. DIVA's correction
    # `u_b = ū/(1 + βF₂)` (Robinson et al. 2022, Eq. 18) is Stage 2 — see
    # `roadmaps/chmy.md`, Phase 3. Reads `depthaverage_x`/`y`, the field the solver
    # actually iterates (decision 1); it used to read `velocity.x`/`y`, which held the same
    # numbers only because `nz == 1` collapsed `ACX3` onto `ACX2`.
    copyto!(asarray(velocity.base_x), asarray(velocity.depthaverage_x))
    copyto!(asarray(velocity.base_y), asarray(velocity.depthaverage_y))
    basalstress!(stress.base_x, stress.base_y, friction.beta_eff,
                velocity.base_x, velocity.base_y, rt, mask)
    return nothing
end

"""
$(TYPEDSIGNATURES)

No-op: `stress.base_x`/`base_y` are untouched, since [`NoFrictionUpdate`](@ref) means
"hold the basal stress fixed at whatever was written there before the solve" — bypassing
the friction law entirely.
"""
update_basalstress!(::MechanicState, ::NoFrictionUpdate, ::Runtime,
                    ::AbstractIceMask = NoMask()) = nothing

"""
$(TYPEDSIGNATURES)

Evaluate the PT velocity rate: velocity gradients, viscosity continuation
([`update_viscosity!`](@ref), a no-op unless `solver.viscosity_continuation` enables it —
must run after the gradients it reads and before the membrane stress that reads its
output, hence its place in this order), membrane stress, basal stress
([`update_basalstress!`](@ref), a no-op if `solver.friction_update` is
[`NoFrictionUpdate`](@ref) — from the SSA-limit basal velocity otherwise; a real DIVA
`F₂` integral is Phase 3 future work), then the PT
velocity rate. Writes `solver.velocity_x_dt`/`velocity_y_dt`, damped by `gamma`, and
`solver.residual_x`/`residual_y`, the undamped rate (see [`dotvel!`](@ref)).

`gamma` is a keyword because [`AutotunedDynamicRelaxation`](@ref) recomputes the damping
*during* the solve, so the value in force at a given iteration is not a property of the
solver; it defaults to `1` (undamped), matching [`dotvel!`](@ref).
"""
# Per-balance grid requirements. SSA is depth-independent by construction: it never reads a
# column field, so any `nz` is fine and a column grid simply costs nothing. DIVA is the
# opposite — it is *defined* by its vertical structure, and on `nz == 1` the single
# quadrature point sits at `σ = ½`, which makes `F₂ = H/(4µ)` against the true `H/(3µ)`
# (see `viscosity_integrals!`). That is 25% low while looking entirely plausible, so it is
# an error rather than a warning (`roadmaps/chmy.md`, Phase 3, decision 4).
_check_momentum_grid(::SSAMomentumBalance, ::Runtime) = nothing

function _check_momentum_grid(::DIVAMomentumBalance, rt::Runtime)
    rt.grid === rt.grid2d && throw(ArgumentError(
        "DIVAMomentumBalance requires a column StaggeredGrid (nz > 1). On a " *
        "depth-averaged grid (nz == 1, rt.grid === rt.grid2d) the viscosity integral F₂ " *
        "is 25% low — a plausible-looking wrong answer, not an approximation. Build the " *
        "grid with a `layering` argument, e.g. " *
        "`StaggeredGrid(T, lx, ly, dx, dy, CorrectedVerticalLayering(T, " *
        "QuadraticSigmaTransform(T, nz)))`, or use SSAMomentumBalance()."))
    return nothing
end

function pseudo_rate!(mech::MechanicState, c::Constants, rt::Runtime,
                      momentum::MomentumBalance2D,
                      solver::PseudoTransientSolver, mask::AbstractIceMask = NoMask();
                      gamma = 1)
    (; velocity, material, topography, stress) = mech

    depthaverage_velocitygradients!(velocity, rt, mask)
    _iterate_viscosity!(mech, momentum, solver, rt, mask)
    membranestress!(stress, velocity, material, topography, momentum, rt, mask)

    update_basalstress!(mech, solver.friction_update, rt, mask)

    dotvel!(solver.velocity_x_dt, solver.velocity_y_dt,
           stress.membrane_xx, stress.membrane_xy, stress.membrane_yy,
           stress.base_x, stress.base_y, stress.driving_x, stress.driving_y,
           topography.thickness, c.density_ice, rt, momentum, mask;
           gamma, resid_x = solver.residual_x, resid_y = solver.residual_y)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Iterate `mech.velocity.x`/`y` in pseudo-time until the SSA/DIVA momentum-balance residual
vanishes. Returns a named tuple
`(; iterations, error, converged, residual, damping, lambda_min)`. `error`/`converged` are
decided by `solver.convergence`; `residual` is the max-norm raw rate `‖r(u)/(ρH)‖` (`solver.residual_x`/`residual_y`, see
[`dotvel!`](@ref)) at the final iterate — always reported, whether or not it also gates the
loop; `damping` is the `γ` in force at the end of the solve and `lambda_min` the last
spectral estimate behind it (`lambda_min` is `NaN` unless `solver.tuning` is
[`AutotunedDynamicRelaxation`](@ref)).

Three solver-held strategies steer the iteration:

 - `solver.pseudo_timestep` ([`AbstractPseudoTimeStep`](@ref)) — how `Δτ` is chosen.
   [`GershgorinPseudoTimeStep`](@ref) (default) bounds the residual operator's spectral
   radius, basal drag included; [`ViscosityPseudoTimeStep`](@ref) is Sandip's `Δτ ∝ 1/η`,
   which **diverges on real geometry with a friction law being solved for** because it
   omits `β/(ρH)`.
 - `solver.convergence` ([`AbstractPTConvergence`](@ref)) — what `abstol` measures.
   [`VelocityIncrement`](@ref) (default) is the exact max-norm increment `|u_new - u_old|`
   (`ux`/`uy` still hold the pre-update iterate in `ux_old`/`uy_old` at that point) rather
   than the algebraic reconstruction `theta_v * dtau * dv` an earlier version used — that
   reconstruction relied on `dtau` being the same scalar everywhere, which a local `Δτ`
   field is not. [`ScaledResidual`](@ref) is the driving-stress-normalized momentum
   residual, and is what a domain mixing drag regimes (grounded ice + ice shelves) needs:
   the increment criterion cannot distinguish a converged shelf from a stalled one.
 - `solver.tuning` ([`AbstractPTTuning`](@ref)) — where `Δτ` and the damping come from.
   [`FixedTuning`](@ref) (default) carries them as hand-set numbers;
   [`AutotunedDynamicRelaxation`](@ref) derives both from a Gershgorin `λ_max` and a
   Rayleigh-quotient `λ_min`, re-derived every `cadence` iterations.

!!! warning "`material.viscosity_depthaveraged`, `topography.thickness` and `friction.beta_eff` need their halo filled by the caller"
    Unlike the velocity, these are fixed inputs for the whole solve, so this function does
    not refresh their halo itself (that is a Phase 4, `bc!`-policy decision belonging to
    whatever step last wrote them, not to the solver reading them). Leaving a halo at its
    allocation default (`setdata!` only touches the interior, by design — see its
    docstring) is not merely inexact at the domain edge, it can be `NaN`: the membrane
    stress's `hlerp(η)` reads the outer ghost ring of `viscosity_depthaveraged` at the two
    boundary `ab` corners, and `hlerp` of an unset (zero) neighbour is `NaN`, exactly like
    the already-documented ice-free case — except here there is no ice-free cell at all,
    only an unfilled halo. Found while writing this solver's own tests, which fill every
    such field with `fill_analytic!` for exactly this reason.

Requires a depth-averaged grid (`rt.grid2d === rt.grid`, i.e. `nz == 1`, checked): DIVA's
vertical-shear integral is Phase 3 future work, so today `mech.velocity.x`/`y` *are* the
depth-averaged velocity being solved for.

Two accelerations over the plain Sandip iteration (`roadmaps/PT-autotune.md` Phase 1):
the pseudo-time step is a *local* field (`solver.dtau_x`/`dtau_y`, via
[`pseudo_dt!`](@ref)) rather than one scalar throttled by the domain's single stiffest
cell, and the velocity rate is damped (via [`dotvel!`](@ref)) rather than recomputed from
scratch every iteration, which is what buys sub-quadratic iteration scaling.
"""
function pseudo_transient!(mech::MechanicState, c::Constants, solver::PseudoTransientSolver,
                           rt::Runtime,
                           momentum::MomentumBalance2D = SSAMomentumBalance(),
                           mask::AbstractIceMask = NoMask())
    _check_momentum_grid(momentum, rt)

    (; velocity) = mech
    (; abstol, maxiter, printout_every, ncheck, tuning) = solver

    ux, uy = velocity.depthaverage_x, velocity.depthaverage_y
    ux_old, uy_old = solver.velocity_x_old, solver.velocity_y_old
    dvx, dvy = solver.velocity_x_dt, solver.velocity_y_dt
    dtau_x, dtau_y = solver.dtau_x, solver.dtau_y
    resid_x, resid_y = solver.residual_x, solver.residual_y

    # Fields held fixed over the PT iteration.
    drivingstress!(mech, c, rt, mask)
    state = _tuning_init!(tuning, solver, mech, c, rt, mask)
    # After `drivingstress!` (it reads the driving stress) and before the loop (it borrows
    # `residual_x`/`residual_y` as scratch, which `dotvel!` overwrites on iteration 1).
    scale = _convergence_scale(solver.convergence, mech, c, solver, rt, mask)

    T = eltype(asarray(ux))
    err  = typemax(T)
    iter = 0
    while err > abstol && iter < maxiter
        iter += 1

        # DIVA's depth-integrated-viscosity chain, at the cadence `solver.div_update` names
        # (never, under the default `NoDIVUpdate`). Before `pseudo_rate!`, so the iteration
        # sees the refreshed `β_eff`/`µ̄` rather than applying them one iteration late.
        #
        # `pseudo_dt!` follows it because `β_eff` feeds the Gershgorin bound behind
        # `dtau_x`/`dtau_y`: a `β_eff` that grew under a stale bound makes that bound
        # optimistic, and the explicit iteration then diverges (decision 8).
        if momentum isa DIVAMomentumBalance && _div_refresh_due(solver.div_update, iter)
            diva_update!(mech, solver, rt, mask)
            pseudo_dt!(solver, mech, c, rt, mask)
        end

        pseudo_rate!(mech, c, rt, momentum, solver, mask; gamma = state.gamma)
        state = _tune!(tuning, state, solver, mech, c, rt, mask, ux, uy)

        # After `_tune!`, which needs the pre-copy `u_old` (see its docstring).
        copyto!(asarray(ux_old), asarray(ux))
        copyto!(asarray(uy_old), asarray(uy))

        pseudo_vel!(asarray(ux), asarray(ux_old), asarray(dvx), asarray(dtau_x),
                    state.theta_v)
        pseudo_vel!(asarray(uy), asarray(uy_old), asarray(dvy), asarray(dtau_y),
                    state.theta_v)

        bc!(rt.arch, rt.grid2d, ux => Neumann())
        bc!(rt.arch, rt.grid2d, uy => Neumann())

        state = _arm_tuning(tuning, state, solver, ux, uy, iter)

        if iter % ncheck == 0 || iter == maxiter
            err = _pt_error(solver.convergence, solver, ux, uy, scale)
        end
        if iter % printout_every == 0
            println("PT iteration $iter: err = $err, gamma = $(state.gamma)")
        end
    end

    # Diagnostic only (does not gate the loop above): the raw, undamped rate norm at the
    # final iterate, independent of `gamma` — see `dotvel!`'s `resid_x`/`resid_y` note.
    residual = max(maximum(abs, asarray(resid_x)), maximum(abs, asarray(resid_y)))

    return (; iterations = iter, error = err, converged = err <= abstol, residual,
              damping = state.gamma, lambda_min = state.lambda_min)
end

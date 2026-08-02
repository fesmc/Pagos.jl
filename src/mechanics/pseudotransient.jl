###############################################################
# Pseudo-transient momentum solver (Sandip et al., 2024)
###############################################################

"""
$(TYPEDSIGNATURES)

Solve the momentum balance of `m` for the velocity field, dispatching on the solver
stored in `m.solver`. For a [`PseudoTransientSolver`](@ref) this runs
[`pseudo_transient!`](@ref).
"""
velocity!(m::Mechanics, c::Constants) = velocity!(m, c, m.solver)

velocity!(m::Mechanics, c::Constants, ::PseudoTransientSolver) = pseudo_transient!(m, c)

"""
$(TYPEDSIGNATURES)

Iterate the velocity field of `m.state` in pseudo-time until the momentum balance
residual vanishes, following the pseudo-transient (PT) method of Sandip et al. (2024).

Each iteration evaluates the PT velocity rate with [`pseudo_rate!`](@ref) — which
dispatches on `m.momentum` (e.g. [`DIVAMomentumBalance`](@ref)) — and relaxes the
velocity with [`pseudo_vel!`](@ref) using the PT time step [`pseudo_dt`](@ref).
Convergence is reached when the max-norm velocity change per iteration drops below
`m.solver.abstol`.

All field updates are KernelAbstractions kernels or broadcasts, so the iteration runs
on whichever backend (CPU/GPU) owns the state arrays.

Returns a named tuple `(; iterations, error, converged)`.
"""
function pseudo_transient!(m::Mechanics, c::Constants)
    (; state, grid, solver) = m
    (; velocity, material, topography, stress) = state
    (; dx, dy) = grid
    (; theta_v, abstol, maxiter, muB, ndim2, dtau_scaling, printout_every) = solver

    (; ncheck) = solver

    ux = view(velocity.x, :, :, 1)
    uy = view(velocity.y, :, :, 1)
    ux_old = solver.velocity_x_old
    uy_old = solver.velocity_y_old
    dvx = solver.velocity_x_dt
    dvy = solver.velocity_y_dt

    # Fields held fixed over the PT iteration
    drivingstress!(stress.driving_x, stress.driving_y, topography.surface,
        topography.thickness, c.density_ice, c.gravity, dx, dy)
    dtau = dtau_scaling *
        pseudo_dt(c.density_ice, dx, dy, material.viscosity_depthaveraged, muB, ndim2)

    err  = typemax(dtau)
    iter = 0
    while err > abstol && iter < maxiter
        iter += 1
        copyto!(ux_old, ux)
        copyto!(uy_old, uy)

        pseudo_rate!(m, c)
        pseudo_vel!(ux, ux_old, dvx, dtau, theta_v)
        pseudo_vel!(uy, uy_old, dvy, dtau, theta_v)

        # All updates above are stream-ordered launches; the host only waits every
        # `ncheck` iterations, when the scalar reduction below drains the stream.
        # u - u_old == theta_v * dtau * dv, so the max-norm velocity change is
        # available from the rate fields without an extra temporary.
        if iter % ncheck == 0 || iter == maxiter
            err = theta_v * dtau * max(maximum(abs, dvx), maximum(abs, dvy))
        end
        if iter % printout_every == 0
            println("PT iteration $iter: err = $err, dtau = $dtau")
        end
    end
    return (; iterations = iter, error = err, converged = err <= abstol)
end

"""
$(TYPEDSIGNATURES)

Evaluate the pseudo-transient velocity rate `m.solver.velocity_x_dt`,
`m.solver.velocity_y_dt` from the current velocity iterate:

 1. velocity gradients ([`velocitygradients!`](@ref)),
 2. scaled strain rate ([`strainrate!`](@ref), dispatching on `m.momentum`),
 3. basal stress ([`basalstress!`](@ref)),
 4. momentum residual divided by the inertial scale ([`dotvel!`](@ref)).

The driving stress `m.state.stress.driving_*` must already be up to date (it does not
change within the PT iteration, so [`pseudo_transient!`](@ref) computes it once).
"""
function pseudo_rate!(m::Mechanics, c::Constants)
    (; state, grid, momentum, solver) = m
    (; velocity, strainrate, material, topography, stress, friction) = state
    (; dx, dy) = grid

    velocitygradients!(velocity, dx, dy)
    strainrate!(strainrate, velocity, material, topography, momentum)

    # TODO: basal velocity from depth-averaged velocity via the F₂ integral (DIVA);
    # for now the basal velocity is taken equal to the depth-averaged one (SSA limit).
    copyto!(velocity.base_x, view(velocity.x, :, :, 1))
    copyto!(velocity.base_y, view(velocity.y, :, :, 1))
    basalstress!(stress.base_x, stress.base_y, friction.beta_eff,
        velocity.base_x, velocity.base_y)

    dotvel!(solver.velocity_x_dt, solver.velocity_y_dt,
        view(strainrate.xx, :, :, 1), view(strainrate.xy, :, :, 1),
        view(strainrate.yy, :, :, 1),
        stress.base_x, stress.base_y, stress.driving_x, stress.driving_y,
        topography.thickness, c.density_ice, dx, dy, momentum)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Compute the rate of pseudo-transient velocity change `dvx`, `dvy` for a
depth-integrated momentum balance (SSA/DIVA):

`` \\partial_\\tau u = \\left( \\nabla \\cdot \\sigma - \\tau_b - \\tau_d \\right) / (\\rho \\, H) ``

where the membrane-stress divergence is formed in-kernel from the scaled strain-rate
components `sxx`, `sxy`, `syy` (already multiplied by the vertically integrated
viscosity, see [`strainrate!`](@ref)). Cells without ice (`H <= 0`) get a zero rate.
"""
function dotvel!(dvx, dvy, sxx, sxy, syy, base_x, base_y, driving_x, driving_y,
    H, density_ice, dx, dy,
    momentum::Union{SSAMomentumBalance, DIVAMomentumBalance})

    backend = get_backend(dvx)
    kernel! = _dotvel_kernel!(backend)
    kernel!(dvx, dvy, sxx, sxy, syy, base_x, base_y, driving_x, driving_y,
        H, density_ice, dx, dy,
        FlatIndexing(1, size(dvx, 1)), FlatIndexing(1, size(dvx, 2));
        ndrange = size(dvx))
    return nothing
end

@kernel function _dotvel_kernel!(dvx, dvy, sxx, sxy, syy, base_x, base_y,
    driving_x, driving_y, H, density_ice, dx, dy, i_idx, j_idx)

    i, j = @index(Global, NTuple)
    im1, ip1, hx = stencil_fd(i, i_idx)
    jm1, jp1, hy = stencil_fd(j, j_idx)
    @inbounds begin
        shear_x = (sxx[ip1, j] - sxx[im1, j]) / (hx * dx) +
                  (sxy[i, jp1] - sxy[i, jm1]) / (hy * dy)
        shear_y = (sxy[ip1, j] - sxy[im1, j]) / (hx * dx) +
                  (syy[i, jp1] - syy[i, jm1]) / (hy * dy)
        Hij = H[i, j]
        if Hij > 0
            dvx[i, j] = (shear_x - base_x[i, j] - driving_x[i, j]) / (density_ice * Hij)
            dvy[i, j] = (shear_y - base_y[i, j] - driving_y[i, j]) / (density_ice * Hij)
        else
            dvx[i, j] = zero(eltype(dvx))
            dvy[i, j] = zero(eltype(dvy))
        end
    end
end

"""
$(TYPEDSIGNATURES)

Pseudo-transient time step based on Sandip et al. (2024), from the ice density `ρ`,
the grid spacings `dx`, `dy`, the depth-averaged viscosity field `mu`, the
bulk-to-shear viscosity ratio `muB` and the numerical-dimensionality constant `ndim`.
"""
function pseudo_dt(ρ, dx, dy, mu, muB, ndim)
    scaling = ρ * dx * dy / (4 * (1 + muB) * ndim)
    # `asarray` is the identity on a plain array, so this is representation-agnostic: the
    # same code path serves the collocated `mu::AbstractArray` and a Chmy-native
    # `mu::AbstractField` without a second method. Routing through it matters on GPU: a
    # bare `maximum(mu::Field)` falls back to generic scalar `getindex` (see `asarray`'s
    # docstring / `roadmaps/chmy.md`).
    return scaling / maximum(asarray(mu))
end

"""
$(TYPEDSIGNATURES)

Relax the velocity field: `v = v_old + theta_v * dotvel * dtau`.
"""
function pseudo_vel!(v, v_old, dotvel, dtau, theta_v)
    @. v = v_old + theta_v * dotvel * dtau
    return nothing
end

###############################################################
# Chmy-native, C-grid staggered pseudo-transient solver
###############################################################
#
# The flagship consumer of the migration: every staggered piece built so far
# (`velocitygradients!`, the membrane-stress `strainrate!`, `basalstress!`,
# `drivingstress!`) is assembled here into the same iterative solve the collocated
# `pseudo_transient!` runs, on `MechanicState` + `Runtime` rather than `Mechanics`
# (`IceSheet`/`Simulation` ownership of a Field-based `Mechanics` bundle is still Phase 2
# future work — see `roadmaps/chmy.md`).
#
# Two things are new at this level, neither of which lives inside a kernel:
#
#  1. **GPU-safe array ops.** `copyto!`, the `v = v_old + θ·dv·dτ` relaxation and the
#     max-norm convergence check are *not* kernels in the collocated code — they are plain
#     `copyto!`/broadcast/`maximum` calls. A `Chmy.Field` is an `AbstractArray` with no
#     `BroadcastStyle`, so those calls compile and give the right answer on CPU via generic
#     scalar `getindex`/`setindex!`, and are wrong or catastrophically slow on GPU (scalar
#     indexing). Every such call below is routed through `asarray` — exactly the case it
#     exists for (see `roadmaps/chmy.md`).
#  2. **`bc!` on the velocity being iterated.** The membrane-stress kernel reads velocity
#     gradients one ring beyond the interior at `ab`, so the velocity's halo must be valid
#     before every `pseudo_rate!` call, not just once at the end — hence a `bc!` refresh
#     each iteration. Phase 4 (mapping `AbstractIndexing` onto real per-equation boundary
#     conditions) has not happened yet, so `Neumann(0)` (zero-gradient) is used as a
#     placeholder at every margin. This is a stand-in, not a physics decision — but it is
#     exact for the uniform-slab solution the tests below check against, since that
#     solution has zero velocity gradient everywhere.

"""
$(TYPEDSIGNATURES)

Chmy-native, local (per-grid-point) pseudo-transient time step (Sandip et al. 2024, Eq. 7):
writes `dtau_x`, `dtau_y` from the depth-averaged viscosity `mu` (`aa`), staggered onto the
velocity's own node class (`acx`/`acy`) by `lerp` — the same arithmetic-mean choice
[`drivingstress!`](@ref) and the staggered [`dotvel!`](@ref) make for the ice thickness `H`.

Distinguished from the collocated, global-scalar [`pseudo_dt`](@ref) by computing a value
*per grid point* rather than reducing over `maximum(mu)`: a single stiff cell no longer
throttles `Δτ` everywhere else in the domain, and no global reduction is needed to compute
it. Only `mu`'s interior is read (`lerp` needs its halo only at the two domain-boundary
`ab`-style corners the membrane stress itself needs — see [`pseudo_transient!`](@ref)'s own
halo warning), so this call has the same halo requirement as the rest of the solve.

`dtau_scaling` folds in the same safety factor `pseudo_transient!` applies to the
collocated, scalar `dtau` (`dtau_scaling * pseudo_dt(...)`), here multiplied in once up
front rather than broadcast over the field afterwards.
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

function pseudo_dt!(solver::PseudoTransientSolver, ::ViscosityPseudoTimeStep,
                    mech::MechanicState, c::Constants, rt::Runtime,
                    ::AbstractIceMask, dx, dy)
    pseudo_dt!(solver.dtau_x, solver.dtau_y, c.density_ice, dx, dy,
               mech.material.viscosity_depthaveraged, solver.muB, solver.ndim2,
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

Fill `solver.dtau_x`/`dtau_y` with `scale / Λ`, where `Λ` is the Gershgorin absolute row sum
of the mass-scaled residual operator (see the source note above). The prefactor is explicit
because the two schemes that use this bound want different ones from the same `Λ`: the
first-order iteration takes `scale = 2·cfl` (`Δτ ≤ 2/λ_max`), while
[`AutotunedDynamicRelaxation`](@ref) takes `scale = Δτ_DR²` — its per-face pseudo-time step
is `Δτ_DR²/Λ` because `Δτ_DR` enters the damped update *twice*, once on the residual and
once on the accumulator (see that type's docstring for the algebra).

`Λ` itself is recoverable from the output as `scale / dtau`, which is how the autotuner
evaluates the Rayleigh quotient's `M`-inner product without storing a third field.
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

Chmy-native, C-grid staggered [`dotvel!`](@ref): the SSA/DIVA pseudo-transient velocity
rate. The membrane-stress divergence is formed by Chmy's `∂x`/`∂y` acting directly on the
already-staggered `sxx` (`aa`), `sxy` (`ab`), `syy` (`aa`) — no interpolation, by the same
node-algebra argument as [`drivingstress!`](@ref) and [`velocitygradients!`](@ref): `∂x`
of `sxx` and `∂y` of `sxy` both land on `acx`; `∂x` of `sxy` and `∂y` of `syy` both land on
`acy`. Only the ice thickness `H` needs staggering onto the face, by `lerp`.

Distinguished from the collocated method by taking a [`Runtime`](@ref) in place of `dx`,
`dy`. Depth-integrated throughout, so it runs on `rt.grid2d`.

The rate is accumulated with damping (Sandip et al. 2024, Eq. 12–14):
`dvx = (1 - gamma) * dvx_old + r_x(u)`, where `dvx_old` is whatever `dvx` already holds on
entry. `gamma = 1` (the default) discards `dvx_old` entirely and reproduces the plain,
undamped rate every call — the collocated method's only behaviour, and this method's
behaviour before `gamma` was added. `0 < gamma < 1` gives the rate memory, which is what
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
                 momentum::Union{SSAMomentumBalance, DIVAMomentumBalance},
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
_convergence_scale(::VelocityIncrement, mech, c, solver, rt, mask) = one(eltype(asarray(mech.velocity.x)))

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
#
# The scheme, the algebra and the warm-up are documented on `AutotunedDynamicRelaxation`
# (`src/mechanics/solvers.jl`); this section is the machinery. Three things are worth
# knowing before reading it.
#
#  1. **Nothing new is stored.** The Rayleigh quotient needs `Δu = u^k - u^{k-1}` and
#     `Δr̃ = r̃(u^k) - r̃(u^{k-1})`, i.e. the residual at two consecutive iterates — and the
#     loop only ever holds one. Rather than allocate a second residual pair for every solver
#     ever constructed, the numerator is split across the two iterations that already have
#     the pieces:
#
#       iteration k,   after `pseudo_vel!`:  u = u^k, u_old = u^{k-1}, resid = r̃(u^{k-1})
#                                            ⟹ Δu·r̃(u^{k-1}) and Δuᵀ M Δu   (`_arm_tuning`)
#       iteration k+1, after `pseudo_rate!`: u, u_old unchanged,       resid = r̃(u^k)
#                                            ⟹ Δu·r̃(u^k)                    (`_tune!`)
#
#     Two scalars carried between them, no fields. This is why `pseudo_transient!` copies
#     `u → u_old` *after* `pseudo_rate!` rather than at the top of the iteration: it is the
#     only reordering that leaves `Δu` intact at both sample points, and it is invisible to
#     everything else (`pseudo_rate!` never reads `u_old`, and `pseudo_vel!` still gets the
#     current iterate as its base).
#
#  2. **`M = diag(Λ)` is recovered from `Δτ`, not stored either.** `gershgorin_dt!` writes
#     `dtau = scale/Λ`, so `Λ = scale/dtau` and `Δuᵀ M Δu = scale · Σ Δu²/dtau`.
#
#  3. **The tuning state is an immutable `NamedTuple`** rebound at those two points, with
#     the same field set for every `AbstractPTTuning`, so the loop is type stable whichever
#     tuning is selected and `FixedTuning`'s methods compile away to nothing.

# `scale`/`lambda_min` are meaningless under `FixedTuning` (nothing reads them); `NaN` says
# so rather than a zero that could be mistaken for a measurement.
_tuning_state(::Type{T}, gamma, theta_v, scale) where {T} =
    (; gamma = T(gamma), theta_v = T(theta_v), scale = T(scale), lambda_min = T(NaN),
       rayleigh_ur = zero(T), rayleigh_uu = zero(T), armed = false)

"""
$(TYPEDSIGNATURES)

Fill `solver.dtau_x`/`dtau_y` for the first iteration and return the initial tuning state
(`gamma`, `theta_v`, the `Δτ` prefactor and the `λ_min` bookkeeping) that
[`pseudo_transient!`](@ref) threads through its loop.

[`FixedTuning`](@ref) defers to [`pseudo_dt!`](@ref) and freezes `solver.gamma`/`theta_v`.
[`AutotunedDynamicRelaxation`](@ref) starts in its warm-up: undamped (`γ = 1`) at the
damped-Jacobi step `Δτ = 1/λ_max`, i.e. `scale = 1` — see that type's docstring for why the
DR step cannot be used before the damping it is paired with is known.
"""
function _tuning_init!(::FixedTuning, solver::PseudoTransientSolver, mech::MechanicState,
                       c::Constants, rt::Runtime, mask::AbstractIceMask)
    pseudo_dt!(solver, mech, c, rt, mask)
    return _tuning_state(typeof(solver.gamma), solver.gamma, solver.theta_v, NaN)
end

function _tuning_init!(::AutotunedDynamicRelaxation, solver::PseudoTransientSolver,
                       mech::MechanicState, c::Constants, rt::Runtime,
                       mask::AbstractIceMask)
    T = typeof(solver.gamma)
    scale = T(solver.dtau_scaling)
    gershgorin_dt!(solver, mech, c, rt, mask, scale)
    return _tuning_state(T, 1, 1, scale)
end

@inline _du_dot_r(u, u_old, r) = (u - u_old) * r

# Guarded because `dtau` is exactly zero off-mask, where `Δu` is zero too: an unguarded
# 0/0 would poison the whole reduction with `NaN`.
@inline _du2_over_dtau(u, u_old, dtau) =
    dtau > zero(dtau) ? (u - u_old)^2 / dtau : zero(dtau)

_sum_du_dot_r(solver, ux, uy) =
    mapreduce(_du_dot_r, +, asarray(ux), asarray(solver.velocity_x_old),
              asarray(solver.residual_x)) +
    mapreduce(_du_dot_r, +, asarray(uy), asarray(solver.velocity_y_old),
              asarray(solver.residual_y))

_sum_du2_over_dtau(solver, ux, uy) =
    mapreduce(_du2_over_dtau, +, asarray(ux), asarray(solver.velocity_x_old),
              asarray(solver.dtau_x)) +
    mapreduce(_du2_over_dtau, +, asarray(uy), asarray(solver.velocity_y_old),
              asarray(solver.dtau_y))

"""
$(TYPEDSIGNATURES)

Sample the half of the Rayleigh quotient that is only available at the end of a PT
iteration — `Δuᵀ r̃(u^{k-1})` and `Δuᵀ M Δu` — and arm [`_tune!`](@ref) to take the other
half at the start of the next one. A no-op under [`FixedTuning`](@ref), and under
[`AutotunedDynamicRelaxation`](@ref) except every `cadence` iterations.

Must be called after `pseudo_vel!`, where `ux`/`uy` hold `u^k` and
`solver.velocity_*_old` hold `u^{k-1}`.
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

Close the Rayleigh quotient armed by [`_arm_tuning`](@ref), derive `Δτ` and the damping
from it (Duretz et al. 2026 Eq. 19–21; the closed form is in
[`AutotunedDynamicRelaxation`](@ref)'s docstring) and refill `solver.dtau_x`/`dtau_y` — the
refill doubling as the `λ_max` re-estimation the nonlinear case needs, since it rebuilds
the Gershgorin bound from the viscosity as it now stands.

Must be called after `pseudo_rate!` and before the `u → u_old` copy, where `ux`/`uy` still
hold `u^k`, `solver.velocity_*_old` still hold `u^{k-1}` and `solver.residual_*` hold
`r̃(u^k)`. A no-op under [`FixedTuning`](@ref), and whenever the quotient is degenerate (a
converged or stationary iterate leaves `Δu = 0`), in which case the previous parameters
stand.
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

    # A Rayleigh quotient of `Â` cannot exceed `λ_max ≤ 1`, so anything above is noise, not
    # a measurement — which is what `Δu` degenerates into once the solve is converged to
    # roundoff and the loop keeps re-estimating from an increment that is pure rounding.
    λ_min = min(T(numerator / denominator), one(T))
    # Duretz Eq. 19 (`c = c_damp·2√λ_min`) and Eq. 20 (`Δτ = c_CFL·2/√λ_max`, with
    # `λ_max ≤ 1` by construction) solved jointly against the damped iteration's exact
    # stability bound `Δτ² = c_CFL²·2(2 - γ)`, so that `γ = c·Δτ` is consistent with the
    # `Δτ` it is paired with rather than with an undamped one.
    # `convert`ed rather than promoted: the tuning and Δτ types are independent keyword
    # arguments of the solver, and a `Float64` `cfl` on a `Float32` solver must not widen
    # the state tuple (which would retype `state` inside the loop).
    d = 2 * convert(T, tu.c_damp) * sqrt(λ_min)
    cc = convert(T, solver.pseudo_timestep.cfl)^2
    dtau = -cc * d + sqrt(cc^2 * d^2 + 4 * cc)
    scale = dtau^2 * convert(T, solver.dtau_scaling)
    gershgorin_dt!(solver, mech, c, rt, mask, scale)
    return merge(state, (; gamma = d * dtau, scale, lambda_min = λ_min, armed = false))
end

###############################################################
# Chmy-native Glen viscosity continuation
###############################################################
#
# `AbstractViscosityContinuation`/`NoViscosityContinuation`/`GlenViscosityContinuation`
# live in `src/mechanics/solvers.jl` next to `PseudoTransientSolver` (the dispatch type they
# extend); `update_viscosity!` lives here because it is PT-loop behaviour, called from
# `pseudo_rate!` below, exactly where the module's other Chmy-native rate machinery lives.

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
                   strainrate.effective, vc.n_glen, vc.strainrate_reg, vc.theta_mu, mask))
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

###############################################################
# Chmy-native basal friction update
###############################################################
#
# `AbstractFrictionUpdate`/`ActiveFrictionUpdate`/`NoFrictionUpdate` live in
# `src/mechanics/solvers.jl` next to `PseudoTransientSolver` (the dispatch type they
# extend), exactly like the viscosity-continuation trio above; `update_basalstress!` lives
# here for the same reason `update_viscosity!` does.

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
    copyto!(asarray(velocity.base_x), asarray(velocity.x))
    copyto!(asarray(velocity.base_y), asarray(velocity.y))
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

Chmy-native [`pseudo_rate!`](@ref): velocity gradients, viscosity continuation
([`update_viscosity!`](@ref), a no-op unless `solver.viscosity_continuation` enables it —
must run after the gradients it reads and before the membrane stress that reads its
output, hence its place in this order), membrane stress, basal stress
([`update_basalstress!`](@ref), a no-op if `solver.friction_update` is
[`NoFrictionUpdate`](@ref) — from the SSA-limit basal velocity otherwise, same `TODO` as
the collocated method: a real DIVA `F₂` integral is Phase 3 future work), then the PT
velocity rate. Writes `solver.velocity_x_dt`/`velocity_y_dt`, damped by `gamma`, and
`solver.residual_x`/`residual_y`, the undamped rate (see [`dotvel!`](@ref)).

`gamma` defaults to `solver.gamma` and is a keyword only because
[`AutotunedDynamicRelaxation`](@ref) recomputes the damping *during* the solve, so the
value in force at a given iteration is not a property of the (immutable) solver.
"""
function pseudo_rate!(mech::MechanicState, c::Constants, rt::Runtime,
                      momentum::Union{SSAMomentumBalance, DIVAMomentumBalance},
                      solver::PseudoTransientSolver, mask::AbstractIceMask = NoMask();
                      gamma = solver.gamma)
    (; velocity, strainrate, material, topography, stress) = mech

    velocitygradients!(velocity, topography.thickness, rt, mask)
    update_viscosity!(mech, solver.viscosity_continuation, rt, mask)
    strainrate!(strainrate, velocity, material, topography, momentum, rt, mask)

    update_basalstress!(mech, solver.friction_update, rt, mask)

    dotvel!(solver.velocity_x_dt, solver.velocity_y_dt,
           strainrate.xx, strainrate.xy, strainrate.yy,
           stress.base_x, stress.base_y, stress.driving_x, stress.driving_y,
           topography.thickness, c.density_ice, rt, momentum, mask;
           gamma, resid_x = solver.residual_x, resid_y = solver.residual_y)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Chmy-native, C-grid staggered [`pseudo_transient!`](@ref): iterate `mech.velocity.x`/`y`
in pseudo-time until the SSA/DIVA momentum-balance residual vanishes — the same algorithm
as the collocated method (see its docstring), with every array-level op routed through
`asarray` and a `Neumann(0)` halo refresh on the velocity each iteration (see the module
note above). Returns a named tuple
`(; iterations, error, converged, residual, damping, lambda_min)`: the first three as the
collocated method does (`error`/`converged` are decided by `solver.convergence`);
`residual` is the max-norm raw rate `‖r(u)/(ρH)‖` (`solver.residual_x`/`residual_y`, see
[`dotvel!`](@ref)) at the final iterate — always reported, whether or not it also gates the
loop; `damping` is the `γ` in force at the end of the solve and `lambda_min` the last
spectral estimate behind it, both `NaN`/`solver.gamma` unless `solver.tuning` is
[`AutotunedDynamicRelaxation`](@ref).

Three solver-held strategies steer the iteration, all defaulting to the behaviour this
method had before they existed, so an unadorned `PseudoTransientSolver(grid)` runs the
original scheme unchanged:

 - `solver.pseudo_timestep` ([`AbstractPseudoTimeStep`](@ref)) — how `Δτ` is chosen.
   [`ViscosityPseudoTimeStep`](@ref) (default) is Sandip's `Δτ ∝ 1/η`;
   [`GershgorinPseudoTimeStep`](@ref) bounds the residual operator's spectral radius
   including the basal-drag term. **On real geometry with a friction law being solved for,
   the default diverges** — Sandip's bound omits `β/(ρH)`, which dominates under grounded
   ice. See [`GershgorinPseudoTimeStep`](@ref).
 - `solver.convergence` ([`AbstractPTConvergence`](@ref)) — what `abstol` measures.
   [`VelocityIncrement`](@ref) (default) is the exact max-norm increment `|u_new - u_old|`
   (`ux`/`uy` still hold the pre-update iterate in `ux_old`/`uy_old` at that point) rather
   than the algebraic reconstruction `theta_v * dtau * dv` an earlier version used — that
   reconstruction relied on `dtau` being the same scalar everywhere, which a local `Δτ`
   field is not. [`ScaledResidual`](@ref) is the driving-stress-normalized momentum
   residual, and is what a domain mixing drag regimes (grounded ice + ice shelves) needs:
   the increment criterion cannot distinguish a converged shelf from a stalled one.
 - `solver.tuning` ([`AbstractPTTuning`](@ref)) — where `Δτ` and the damping come from.
   [`FixedTuning`](@ref) (default) uses the hand-set `solver.theta_v`/`solver.gamma`;
   [`AutotunedDynamicRelaxation`](@ref) derives both from a Gershgorin `λ_max` and a
   Rayleigh-quotient `λ_min`, re-derived every `cadence` iterations, and ignores those two
   fields.

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
depth-averaged velocity being solved for, exactly as in the collocated method (which reads
`view(velocity.x, :, :, 1)`).

Two accelerations over the collocated method, both from Sandip et al. (2024) and both
solver-level parameters (`solver.gamma`, `roadmaps/PT-autotune.md` Phase 1):

 - the pseudo-time step is a *local* field (`solver.dtau_x`/`dtau_y`, via
   [`pseudo_dt!`](@ref)) rather than one scalar throttled by the domain's single stiffest
   cell;
 - the velocity rate is damped (`solver.gamma`, via [`dotvel!`](@ref)) rather than fully
   recomputed from scratch every iteration, which is what buys sub-quadratic iteration
   scaling. `gamma = 1` (the default) disables damping and recovers the plain iteration.
"""
function pseudo_transient!(mech::MechanicState, c::Constants, solver::PseudoTransientSolver,
                           rt::Runtime,
                           momentum::Union{SSAMomentumBalance, DIVAMomentumBalance} = DIVAMomentumBalance(),
                           mask::AbstractIceMask = NoMask())
    rt.grid === rt.grid2d || throw(ArgumentError(
        "the Chmy-native pseudo_transient! only supports a depth-averaged StaggeredGrid " *
        "(nz == 1, i.e. rt.grid === rt.grid2d); DIVA's vertical shear integral is not " *
        "yet ported (roadmaps/chmy.md, Phase 3)."))

    (; velocity) = mech
    (; abstol, maxiter, printout_every, ncheck, tuning) = solver

    ux, uy = velocity.x, velocity.y
    ux_old, uy_old = solver.velocity_x_old, solver.velocity_y_old
    dvx, dvy = solver.velocity_x_dt, solver.velocity_y_dt
    dtau_x, dtau_y = solver.dtau_x, solver.dtau_y
    resid_x, resid_y = solver.residual_x, solver.residual_y

    # Fields held fixed over the PT iteration.
    drivingstress!(mech, c, rt, mask)
    # Fills `dtau_x`/`dtau_y` and hands back the iteration parameters (see the autotuning
    # section above); under the default `FixedTuning` this is `pseudo_dt!` plus the
    # solver's own frozen `gamma`/`theta_v`.
    state = _tuning_init!(tuning, solver, mech, c, rt, mask)
    # After `drivingstress!` (it reads the driving stress) and before the loop (it borrows
    # `residual_x`/`residual_y` as scratch, which `dotvel!` overwrites on iteration 1).
    scale = _convergence_scale(solver.convergence, mech, c, solver, rt, mask)

    T = eltype(asarray(ux))
    err  = typemax(T)
    iter = 0
    while err > abstol && iter < maxiter
        iter += 1

        pseudo_rate!(mech, c, rt, momentum, solver, mask; gamma = state.gamma)
        # Between the rate and the copy on purpose: this is the only point where the
        # residual at the current iterate and the *previous* increment coexist, which is
        # what lets the autotuner close its Rayleigh quotient without a second residual
        # buffer (see the autotuning section above). A no-op under `FixedTuning`.
        state = _tune!(tuning, state, solver, mech, c, rt, mask, ux, uy)

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

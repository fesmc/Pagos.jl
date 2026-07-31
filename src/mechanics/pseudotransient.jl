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

"""
$(TYPEDSIGNATURES)

Chmy-native [`pseudo_rate!`](@ref): velocity gradients, viscosity continuation
([`update_viscosity!`](@ref), a no-op unless `solver.viscosity_continuation` enables it —
must run after the gradients it reads and before the membrane stress that reads its
output, hence its place in this order), membrane stress, basal stress (from the SSA-limit
basal velocity — same `TODO` as the collocated method: a real DIVA `F₂` integral is Phase 3
future work), then the PT velocity rate. Writes `solver.velocity_x_dt`/`velocity_y_dt`,
damped by `solver.gamma`, and `solver.residual_x`/`residual_y`, the undamped rate (see
[`dotvel!`](@ref)).
"""
function pseudo_rate!(mech::MechanicState, c::Constants, rt::Runtime,
                      momentum::Union{SSAMomentumBalance, DIVAMomentumBalance},
                      solver::PseudoTransientSolver, mask::AbstractIceMask = NoMask())
    (; velocity, strainrate, material, topography, stress, friction) = mech

    velocitygradients!(velocity, topography.thickness, rt, mask)
    update_viscosity!(mech, solver.viscosity_continuation, rt, mask)
    strainrate!(strainrate, velocity, material, topography, momentum, rt, mask)

    copyto!(asarray(velocity.base_x), asarray(velocity.x))
    copyto!(asarray(velocity.base_y), asarray(velocity.y))
    basalstress!(stress.base_x, stress.base_y, friction.beta_eff,
                velocity.base_x, velocity.base_y, rt, mask)

    dotvel!(solver.velocity_x_dt, solver.velocity_y_dt,
           strainrate.xx, strainrate.xy, strainrate.yy,
           stress.base_x, stress.base_y, stress.driving_x, stress.driving_y,
           topography.thickness, c.density_ice, rt, momentum, mask;
           gamma = solver.gamma, resid_x = solver.residual_x, resid_y = solver.residual_y)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Chmy-native, C-grid staggered [`pseudo_transient!`](@ref): iterate `mech.velocity.x`/`y`
in pseudo-time until the SSA/DIVA momentum-balance residual vanishes — the same algorithm
as the collocated method (see its docstring), with every array-level op routed through
`asarray` and a `Neumann(0)` halo refresh on the velocity each iteration (see the module
note above). Returns a named tuple `(; iterations, error, converged, residual)`: the first
three as the collocated method does (`error`/`converged` are decided from the max-norm
velocity increment); `residual` is the additional max-norm raw rate `‖r(u)/(ρH)‖`
(`solver.residual_x`/`residual_y`, see [`dotvel!`](@ref)) at the final iterate — a
diagnostic only, not part of the stopping condition (see `roadmaps/PT-autotune.md`, Phase 1,
for why: unifying the two would change what `abstol` means, which needs its own pass).

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

The convergence check is the exact max-norm velocity increment `|u_new - u_old|`
(`ux`/`uy` still hold the pre-update iterate in `ux_old`/`uy_old` at that point) rather
than an algebraic reconstruction from `theta_v * dtau * dv` — that reconstruction relied on
`dtau` being the same scalar everywhere, which a local `Δτ` field no longer is.
"""
function pseudo_transient!(mech::MechanicState, c::Constants, solver::PseudoTransientSolver,
                           rt::Runtime,
                           momentum::Union{SSAMomentumBalance, DIVAMomentumBalance} = DIVAMomentumBalance(),
                           mask::AbstractIceMask = NoMask())
    rt.grid === rt.grid2d || throw(ArgumentError(
        "the Chmy-native pseudo_transient! only supports a depth-averaged StaggeredGrid " *
        "(nz == 1, i.e. rt.grid === rt.grid2d); DIVA's vertical shear integral is not " *
        "yet ported (roadmaps/chmy.md, Phase 3)."))

    (; velocity, material) = mech
    (; theta_v, abstol, maxiter, muB, ndim2, dtau_scaling, printout_every, ncheck) = solver

    ux, uy = velocity.x, velocity.y
    ux_old, uy_old = solver.velocity_x_old, solver.velocity_y_old
    dvx, dvy = solver.velocity_x_dt, solver.velocity_y_dt
    dtau_x, dtau_y = solver.dtau_x, solver.dtau_y
    resid_x, resid_y = solver.residual_x, solver.residual_y

    dx = Δx(rt.grid2d, Center(), 1, 1, 1)
    dy = Δy(rt.grid2d, Center(), 1, 1, 1)

    # Fields held fixed over the PT iteration.
    drivingstress!(mech, c, rt, mask)
    pseudo_dt!(dtau_x, dtau_y, c.density_ice, dx, dy, material.viscosity_depthaveraged,
              muB, ndim2, dtau_scaling, rt)

    T = eltype(asarray(ux))
    err  = typemax(T)
    iter = 0
    abs_diff(a, b) = abs(a - b)
    while err > abstol && iter < maxiter
        iter += 1
        copyto!(asarray(ux_old), asarray(ux))
        copyto!(asarray(uy_old), asarray(uy))

        pseudo_rate!(mech, c, rt, momentum, solver, mask)
        pseudo_vel!(asarray(ux), asarray(ux_old), asarray(dvx), asarray(dtau_x), theta_v)
        pseudo_vel!(asarray(uy), asarray(uy_old), asarray(dvy), asarray(dtau_y), theta_v)

        bc!(rt.arch, rt.grid2d, ux => Neumann())
        bc!(rt.arch, rt.grid2d, uy => Neumann())

        if iter % ncheck == 0 || iter == maxiter
            err = max(mapreduce(abs_diff, max, asarray(ux), asarray(ux_old)),
                      mapreduce(abs_diff, max, asarray(uy), asarray(uy_old)))
        end
        if iter % printout_every == 0
            println("PT iteration $iter: err = $err")
        end
    end

    # Diagnostic only (does not gate the loop above): the raw, undamped rate norm at the
    # final iterate, independent of `gamma` — see `dotvel!`'s `resid_x`/`resid_y` note.
    residual = max(maximum(abs, asarray(resid_x)), maximum(abs, asarray(resid_y)))

    return (; iterations = iter, error = err, converged = err <= abstol, residual)
end

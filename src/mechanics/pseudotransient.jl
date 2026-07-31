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
    return scaling / maximum(mu)
end

"""
$(TYPEDSIGNATURES)

Relax the velocity field: `v = v_old + theta_v * dotvel * dtau`.
"""
function pseudo_vel!(v, v_old, dotvel, dtau, theta_v)
    @. v = v_old + theta_v * dotvel * dtau
    return nothing
end

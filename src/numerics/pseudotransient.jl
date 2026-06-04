function update_dynamics!(
    states::Vector{State{T}},
    grids::Vector{RegularGrid{T}},
    dynamics::Vector{<:AbstractDynamics{T}},
    dyn_solvers::Vector{<:AbstractDynamicsSolver{T}},
    tools::Vector{Tools{T}},
    constants::Vector{PhysicalConstants{T}},
    options::Vector{Options},
) where {T<:AbstractFloat}
    for i in 1:ngrids
        
        if dyn_solvers[i] isa LinearSolver
            throw(error("`update_dynamics` not implemented for `LinearSolver` yet"))
        end

        if i > 1
            interpolate!(fine, coarse)
        end

        update_dynamics!(states[i], grids[i], dynamics[i], dyn_solver, tools[i],
            constants[i], options[i])

        interpolate!(common, grids[i])
    end

    return nothing
end

function update_dynamics!(
    state::State{T},
    grid::RegularGrid{T},
    dynamics::DIVA{T},
    dyn_solver::PseudoTransientSolver{T},
    tools::Tools{T},
    constants::PhysicalConstants{T},
    options::Options,
) where {T<:AbstractFloat}
    (; abstol, maxiter, err) = dyn_solver
    (; buffer2D) = tools

    # Init PT loop
    @. err = Inf
    iter = 0

    while (err[max(iter, 1)] > abstol) && (iter + 1 < maxiter)
        iter += 1
        pseudo_update_dynamics!(state, grid, dynamics, dyn_solver, tools, constants,
            options)
        pseudo_transient_error!(err, buffer2D, iter, state.dynamics, dynamics,
            state.masks.update_ice)
        pseudo_transient_monitor(iter, dyn_solver, err, state.dynamics)
    end
    return nothing
end

"""
    pseudo_transient_monitor(iter, dyn_solver, err, state_dynamics)

Print pseudo-transient iteration information.
"""
function pseudo_transient_monitor(iter, dyn_solver, err, state_dynamics)
    if iter % dyn_solver.printout_every == 0
        println("Pseudo-transient iteration: $iter")
        println("Pseudo-transient error: $(err[iter])")
        @show extrema(state_dynamics.v_x_dt)
        @show extrema(state_dynamics.v_old_x)
        @show extrema(state_dynamics.v_x)
        println("------------------")
    end
end

"""
    pseudo_transient_error(state_dynamics, dynamics<:AbstractDynamics)

Calculate the pseudo-transient error depending on the dimensionality of the dynamics.
"""
function pseudo_transient_error!(
    err,
    buffer2D,
    iter,
    state_dynamics::DynamicsState2D{T},
    dynamics::Dynamics2D{T},
    mask::B,
) where {B}

    buffer2D .= 0
    @inbounds for I in view(mask)
        buffer2D[I] = abs(state_dynamics.v_x[I] - state_dynamics.v_old_x[I]) +
            abs(state_dynamics.v_y[I] - state_dynamics.v_old_y[I])
    end
    err[iter] = maximum(buffer2D)

    # buffer2D .= abs.(state_dynamics.v_x .- state_dynamics.v_old_x)
    # err[iter] = maximum(buffer2D)
    # buffer2D .= abs.(state_dynamics.v_y .- state_dynamics.v_old_y)
    # err[iter] = max(err[iter], maximum(buffer2D))

    # @show err[iter]
    return nothing
end

function pseudo_transient_error(state_dynamics, dynamics::D) where
    {T<:AbstractFloat, D<:Dynamics3D{T}}
    return max(
        maximum(abs.(state_dynamics.v_x - state_dynamics.v_old_x)),
        maximum(abs.(state_dynamics.v_y - state_dynamics.v_old_y)),
        maximum(abs.(state_dynamics.v_z - state_dynamics.v_old_z)),
    )
end

"""
    pseudo_update_dynamics!(icesheet::IceSheet)

Perform a pseudo time step to update the velocity field.
"""
function pseudo_update_dynamics!(
    state::State{T},
    grid::RegularGrid{T},
    dyn::AbstractDynamics{T},
    dyn_solver::PseudoTransientSolver{T},
    tools::Tools{T},
    constants::PhysicalConstants{T},
    options::Options,
) where {T<:AbstractFloat}

    # Unpack structs
    (; v_x, v_y, v_old_x, v_old_y, v_x_dt, v_y_dt, dtau) = state.dynamics
    (; theta_v) = dyn_solver

    # PT loop following Sandip et al. (2024)
    v_old_x .= v_x
    v_old_y .= v_y

    # Pre-updating dynamics dispatches on Dynamics2D and Dynamics3D
    pre_update_dynamics!(state, grid, dyn, dyn_solver, tools, constants, options)
    pseudo_vel!(v_x, v_old_x, v_x_dt, dtau, theta_v, state.masks.update_ice)
    pseudo_vel!(v_y, v_old_y, v_y_dt, dtau, theta_v, state.masks.update_ice)

    return nothing
end

function pre_update_dynamics!(
    state::State{T},
    grid::RegularGrid{T},
    dyn::Dynamics2D{T},
    dyn_solver::PseudoTransientSolver{T},
    tools::Tools{T},
    constants::PhysicalConstants{T},
    options::Options,
) where {T<:AbstractFloat}

    # Unpack structs
    (; dynamics, topography, material, masks) = state
    (; v_x, v_y, v_x_dt, v_y_dt, v_x_dx, v_x_dy, v_y_dx, v_y_dy) = dynamics
    (; strain_dt_d2x, strain_dt_dxdy, strain_dt_d2y) = dynamics
    (; stress_shear_x, stress_shear_y, v_basal_x, v_basal_y) = dynamics
    (; stress_basal_x, stress_basal_y, β, β_ac_x, β_ac_y) = dynamics
    (; stress_driving_x, stress_driving_y, dtau) = dynamics

    (; H, z_bed) = topography
    (; mu, N_ab) = material

    (; dx, dy, nx, ny) = grid
    (; dtau_scaling, μ_B, ndim2) = dyn_solver
    (; ρ_ice, g) = constants
    (; buffer2D) = tools

    map_aa2ac!(β_ac_x, β_ac_y, β, masks.ice)
    # check_nans(options, [β_ac_x, β_ac_y], ["β_ac_y, β_ac_x"])

    vintegrated_viscosity!(N_ab, mu, H, nx, ny, masks.ice, masks.ice)
    # check_nans(options, N_ab, "N_ab")

    velocitygradients!(v_x_dx, v_x_dy, v_y_dx, v_y_dy, v_x, v_y, dx, dy, masks.update_ice,
        masks.update_ice)
    # check_nans(options, [v_x_dx, v_x_dy, v_y_dx, v_y_dy], ["v_x_dx", "v_x_dy", "v_y_dx",
    #     "v_y_dy"])

    scaledstrainrate!(strain_dt_d2x, strain_dt_dxdy, strain_dt_d2y, v_x_dx, v_x_dy,
        v_y_dx, v_y_dy, N_ab, masks.update_ice)
    # check_nans(options, [strain_dt_d2x, strain_dt_dxdy, strain_dt_d2y], ["strain_dt_d2x", 
    #     "strain_dt_dxdy", "strain_dt_d2y"])

    stress_shear!(stress_shear_x, stress_shear_y, strain_dt_d2x, strain_dt_dxdy,
        strain_dt_d2y, buffer2D, dx, dy, masks.update_ice, masks.update_ice)
    # check_nans(options, [stress_shear_x, stress_shear_y], ["stress_shear_x", "stress_shear_y"])

    # TODO: basalvelocity!(state, domain, params, options)
    # basal_velocity_from_depthavg_velocity!(v_basal_x, v_x, β, F_2, mask)
    # basal_velocity_from_depthavg_velocity!(v_basal_y, v_y, β, F_2, mask)
    v_basal_x .= v_x
    v_basal_y .= v_y

    stress_basal!(stress_basal_x, stress_basal_y, β_ac_x, β_ac_y, v_basal_x, v_basal_y,
        masks.update_ice)
    # check_nans(options, [stress_basal_x, stress_basal_y], ["stress_basal_x",
    #    "stress_basal_y"])

    stress_driving!(stress_driving_x, stress_driving_y, buffer2D, ρ_ice, g, H, z_bed, dx,
        dy, masks.update_ice, masks.update_ice)
    # check_nans(options, [stress_driving_x, stress_driving_y], ["stress_driving_x",
    #     "stress_driving_y"])

    dotvel!(v_x_dt, stress_shear_x, stress_basal_x, stress_driving_x, ρ_ice, H,
        masks.update_ice, dyn)
    dotvel!(v_y_dt, stress_shear_y, stress_basal_y, stress_driving_y, ρ_ice, H,
        masks.update_ice, dyn)
    # check_nans(options, [v_x_dt, v_y_dt], ["v_x_dt", "v_y_dt"])

    @. dtau = dtau_scaling * pseudo_dt(ρ_ice, dx, dy, mu, μ_B, ndim2)
    # check_nans(options, dtau, "dtau")
    return nothing
end

# 3D case
function pre_update_dynamics!(
    state::State{T},
    grid::RegularGrid{T},
    dyn::D,
    dyn_solver::PseudoTransientSolver{T},
    tools::Tools{T},
    constants::PhysicalConstants{T},
    options::Options,
) where {T<:AbstractFloat, D<:Dynamics3D{T}}
    return nothing
end

function vintegrated_viscosity!(N_ab, mu, H, nx, ny, mask_defined, mask_update)
    @inbounds for I in view(mask_update)
        ip1 = periodic_plusindex(I.I[1], nx)
        jp1 = periodic_plusindex(I.I[2], ny)
        N_ab[I] = 0.25 * (
            H[I]*mu[I] +
            H[ip1, I.I[2]] * mu[ip1, I.I[2]] +
            H[I.I[1], jp1] * mu[I.I[1], jp1] +
            H[ip1, jp1] * mu[ip1, jp1]
        )
    end
end

"""
    dotvel!(dotvel, stress_shear, stress_basal, stress_driving, ρ_ice, H, nx, ny)

Calculate the rate of pseudo-transient velocity change.
"""
function dotvel!(dotvel, stress_shear, stress_basal, stress_driving, ρ_ice,
    H, mask, dyn::DIVA{T}) where {T<:AbstractFloat}
    for I in view(mask)
        dotvel[I] = (stress_shear[I] - stress_basal[I] - stress_driving[I]) /
                (ρ_ice * (H[I] + 1e-3))
    end
    return nothing
end

"""
    pseudo_dt(ρ, dx, mu, μ_B)

Calculate the pseudo-transient time step based on Sandip et al. (2024).
"""
function pseudo_dt(ρ, dx, dy, mu, μ_B, ndim)
    scaling = ρ * dx * dy / (4 * (1 + μ_B) * ndim)
    return minimum( scaling ./ mu )
end

"""
    pseudo_vel!(v, v_old, pseudo_dotvel, dtau, theta_v)

Update the pseudo-transient velocity field.
"""
function pseudo_vel!(v, v_old, pseudo_dotvel, dtau, theta_v)
    @. v = pseudo_vel.(v_old, pseudo_dotvel, dtau, theta_v)
    return nothing
end

function pseudo_vel!(v, v_old, pseudo_dotvel, dtau, theta_v, mask)
    @inbounds for I in view(mask)
        v[I] = pseudo_vel(v_old[I], pseudo_dotvel[I], dtau[I], theta_v)
    end
    return nothing
end

pseudo_vel(v_old, pseudo_dotvel, dtau, theta_v) = v_old + theta_v * pseudo_dotvel * dtau
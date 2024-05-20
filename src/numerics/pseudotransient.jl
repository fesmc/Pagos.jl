"""
    pseudo_dotvel!(icesheet::IceSheet)
    pseudo_dotvel!(state::State, domain::Domain, params::Params, options::Options)

Calculate the pseudo-transient velocity field. This is an implementation
of the pseudo-transient method, where the velocity field is updated until the
change in the velocity field is below a certain threshold. The method is based
on the work of Sandip et al. (2024).
"""
function pseudo_dotvel!(icesheet::IceSheet)
    (; state, domain, params, options) = icesheet
    pseudo_dotvel!(state, domain, params, options)
    return nothing
end

function pseudo_dotvel!(state::State{T}, domain::Domain{T}, params::Params{T},
    options::Options{T}) where {T<:AbstractFloat}
    (; ux, uy, ux_b, uy_b) = state
    (; ux_x, ux_y, uy_x, uy_y) = state
    (; H, z_b, mu, N_ab) = state
    (; strainrate_xx, strainrate_xy, strainrate_yy) = state
    (; shearstress_x, shearstress_y, dotvel_x, dotvel_y) = state
    (; basalstress_x, basalstress_y, beta_acx, beta_acy) = state
    (; drivingstress_x, drivingstress_y, prealloc) = state
    (; dx, dy, nx, ny) = domain
    (; rho_ice, g) = params

    stagger_beta!(state, domain)
    vintegrated_viscosity!(N_ab, prealloc, mu, H, nx, ny)

    velocitygradients!(ux_x, ux_y, uy_x, uy_y, ux, uy, dx, dy, nx, ny)
    if options.debug
        if hasnan(ux_x) || hasnan(ux_y) || hasnan(uy_x) || hasnan(uy_y)
            throw(ArgumentError("NaNs in velocity gradients"))
        end
    end

    scaledstrainrate!(strainrate_xx, strainrate_xy, strainrate_yy, ux_x, ux_y,
        uy_x, uy_y, N_ab)
    if options.debug
        if hasnan(strainrate_xx) || hasnan(strainrate_xy) || hasnan(strainrate_yy)
            throw(ArgumentError("NaNs in strain rate"))
        end
    end

    shearstress!(shearstress_x, shearstress_y, strainrate_xx, strainrate_xy, strainrate_yy,
        prealloc, dx, dy, nx, ny)
    if options.debug
        if hasnan(shearstress_x) || hasnan(shearstress_y)
            throw(ArgumentError("NaNs in shear stress"))
        end
    end

    # TODO: basalvelocity!(state, domain, params, options)
    ux_b .= ux
    uy_b .= uy
    basalstress!(basalstress_x, basalstress_y, beta_acx, beta_acy, ux_b, uy_b)
    if options.debug
        if hasnan(basalstress_x) || hasnan(basalstress_y)
            throw(ArgumentError("NaNs in basal stress"))
        end
    end

    drivingstress!(drivingstress_x, drivingstress_y, prealloc, rho_ice, g, H, z_b, dx, dy, nx, ny)
    if options.debug
        if hasnan(drivingstress_x) || hasnan(drivingstress_y)
            throw(ArgumentError("NaNs in driving stress"))
        end
    end

    dotvel!(dotvel_x, shearstress_x, basalstress_x, drivingstress_x, rho_ice, H, nx, ny)
    dotvel!(dotvel_y, shearstress_y, basalstress_y, drivingstress_y, rho_ice, H, nx, ny)
    if options.debug
        if hasnan(dotvel_x) || hasnan(dotvel_y)
            throw(ArgumentError("NaNs in dotvel"))
        end
    end

    return nothing
end

function vintegrated_viscosity!(N_ab, prealloc, mu, H, nx, ny)
    @. prealloc = H * mu
    for i in 1:nx, j in 1:ny
        ip1 = periodic_bc_plusindex(i, nx)
        jp1 = periodic_bc_plusindex(j, ny)
        N_ab[i, j] = 0.25 * (prealloc[i, j] + prealloc[ip1, j] + prealloc[i, jp1] +
            prealloc[ip1, jp1])
    end
end

"""
    dotvel!(dotvel, shearstress, basalstress, drivingstress, rho_ice, H, nx, ny)

Calculate the rate of pseudo-transient velocity change.
"""
function dotvel!(dotvel, shearstress, basalstress, drivingstress, rho_ice::T,
    H, nx, ny) where {T<:AbstractFloat}

    @inbounds for i in 1:nx, j in 1:ny
        if H[i, j] > 0
            # dotvel[i, j] = ( shearstress[i, j] - basalstress[i, j] ) /
            #     (rho_ice * H[i, j]) - drivingstress[i, j]
            dotvel[i, j] = (shearstress[i, j] - basalstress[i, j] - drivingstress[i, j]) /
                (rho_ice * H[i, j])
        else
            dotvel[i, j] = T(0.0)
        end
    end
    return nothing
end

"""
    pseudo_dt(rho, dx, mu, muB)

Calculate the pseudo-transient time step based on Sandip et al. (2024).
"""
function pseudo_dt(rho, dx, mu, muB)
    scaling = rho * dx^2 / (4 * (1 + muB) * 4.1)
    return minimum( scaling ./ mu )
end

"""
    pseudo_vel!(v, v_old, pseudo_dotvel, dtau, theta_v)

Update the pseudo-transient velocity field.
"""
function pseudo_vel!(v, v_old, pseudo_dotvel, dtau, theta_v)
    @. v = v_old + theta_v * pseudo_dotvel * dtau
    return nothing
end

"""
    pseudo_transient!(icesheet::IceSheet)

Perform the pseudo-transient method to update the velocity field.
"""
function pseudo_transient!(icesheet::IceSheet{T}) where {T<:AbstractFloat}
    # Unpack structs
    (; state, domain, params, options) = icesheet
    (; ux, uy, ux_old, uy_old, ux, uy) = state
    (; dotvel_x, dotvel_y, mu) = state
    (; rho_ice, muB) = params
    (; dx, dy, nx, ny) = domain

    (; theta_v, abstol, maxiter, dtau_scaling) = options

    # Init PT loop
    err = fill(2 * abstol, maxiter)
    iter = 0
    dtau = 1e-5

    # PT loop following Sandip et al. (2024)
    while err[max(iter, 1)] > abstol && iter + 1 < maxiter
        iter += 1
        ux_old .= ux
        uy_old .= uy

        pseudo_dotvel!(state, domain, params, options)

        dtau = dtau_scaling * pseudo_dt(rho_ice, dx, mu, muB)

        pseudo_vel!(ux, ux_old, dotvel_x, dtau, theta_v)
        pseudo_vel!(uy, uy_old, dotvel_y, dtau, theta_v)

        
        err[iter] = max(
            maximum(abs.(ux - ux_old)),
            maximum(abs.(uy - uy_old)),
        )

        if iter % options.printout_every == 0
            @show extrema(dotvel_x)
            @show extrema(ux_old)
            @show extrema(ux)

            @show iter
            @show dtau
            @show err[iter]
            println("------------------")
        end
    end
    
    return nothing
end
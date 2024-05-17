function pseudo_dotvel!(icesheet::IceSheet)
    (; state, domain, params, options) = icesheet
    pseudo_dotvel!(state, domain, params, options)
    return nothing
end

function pseudo_dotvel!(state::State{T}, domain::Domain{T}, params::Params{T},
    options::Options{T}) where {T<:AbstractFloat}
    (; ux, uy, ux_x, ux_y, uy_x, uy_y, H, z_b, mu) = state
    (; strainrate_xx, strainrate_xy, strainrate_yy) = state
    (; shearstress_x, shearstress_y, dotvel_x, dotvel_y) = state
    (; basalstress_x, basalstress_y, β_acx, β_acy) = state
    (; drivingstress_x, drivingstress_y, prealloc) = state
    (; dx, dy, nx, ny) = domain
    (; rho_ice, g) = params

    velocitygradients!(ux_x, ux_y, uy_x, uy_y, ux, uy, dx, dy, nx, ny)
    if options.debug
        if hasnan(ux_x) || hasnan(ux_y) || hasnan(uy_x) || hasnan(uy_y)
            throw(ArgumentError("NaNs in velocity gradients"))
        end
    end

    strainrate!(strainrate_xx, strainrate_xy, strainrate_yy, ux_x, ux_y, uy_x, uy_y, mu, H)
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

    basalstress!(basalstress_x, basalstress_y, β_acx, β_acy, ux, uy)
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

    # @. dotvel_x = ( shearstress_x - basalstress_x ) / (rho_ice * H) - drivingstress_x
    # @. dotvel_y = ( shearstress_y - basalstress_y ) / (rho_ice * H) - drivingstress_y

    dotvel!(dotvel_x, shearstress_x, basalstress_x, drivingstress_x, rho_ice, H, nx, ny)
    dotvel!(dotvel_y, shearstress_y, basalstress_y, drivingstress_y, rho_ice, H, nx, ny)
    if options.debug
        if hasnan(dotvel_x) || hasnan(dotvel_y)
        throw(ArgumentError("NaNs in dotvel"))
        end
    end

    return nothing
end

function dotvel!(dotvel, shearstress, basalstress, drivingstress, rho_ice::T,
    H, nx, ny) where {T<:AbstractFloat}

    @inbounds for i in 1:nx, j in 1:ny
        if H[i, j] > 0
            dotvel[i, j] = ( shearstress[i, j] - basalstress[i, j] ) /
                (rho_ice * H[i, j]) - drivingstress[i, j]
        else
            dotvel[i, j] = T(0.0)
        end
    end
    return nothing
end

function pseudo_dt(rho, dx, mu, muB)
    scaling = rho * dx^2 / (4 * (1 + muB) * 4.1)
    return minimum( scaling ./ mu )
end

function pseudo_vel!(v, v_old, pseudo_dotvel, dtau, theta_v)
    @. v = v_old + theta_v * pseudo_dotvel * dtau
    return nothing
end

function pseudo_transient!(icesheet::IceSheet{T}) where {T<:AbstractFloat}
    # Unpack structs
    (; state, domain, params, options) = icesheet
    (; pseudo_ux, pseudo_uy, pseudo_ux_old, pseudo_uy_old, ux, uy) = state
    (; dotvel_x, dotvel_y, mu) = state
    (; rho_ice, muB) = params
    (; dx, dy, nx, ny) = domain

    (; theta_v, abstol, maxiter) = options

    # Assign initial values
    pseudo_ux .= ux
    pseudo_uy .= uy

    # Init PT loop
    err = 2 * abstol
    iter = 0
    dtau = 1e-5

    # PT loop following Sandip et al. (2024)
    while err > abstol && iter + 1 < maxiter
        iter += 1
        pseudo_ux_old .= pseudo_ux
        pseudo_uy_old .= pseudo_uy

        pseudo_dotvel!(state, domain, params, options)
        dtau = pseudo_dt(rho_ice, dx, mu, muB)

        pseudo_vel!(pseudo_ux, pseudo_ux_old, dotvel_x, dtau, theta_v)
        pseudo_vel!(pseudo_uy, pseudo_uy_old, dotvel_y, dtau, theta_v)

        if iter % options.compute_residual_every == 0
            err = (norm((pseudo_ux - pseudo_ux_old)) + norm((pseudo_uy - pseudo_uy_old))) / (nx * ny)
            @show iter
            @show dtau
            @show err
        end

    end
    
    return nothing
end
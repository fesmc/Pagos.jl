function pseudo_dotvel!(state, domain, params)
    (; ux, uy, ux_x, ux_y, uy_x, uy_y, H, z_b, mu) = state
    (; strainrate_xx, strainrate_xy, strainrate_yy) = state
    (; shearstress_x, shearstress_y, pseudo_ux, pseudo_uy) = state
    (; basalstress_x, basalstress_y, β_acx, β_acy) = state
    (; drivingstress_x, drivingstress_y) = state
    (; dx, dy, nx, ny) = domain
    (; rho_ice, g) = params
    velocitygradients!(ux_x, ux_y, uy_x, uy_y, ux, uy, dx, dy, nx, ny)
    strainrate!(strainrate_xx, strainrate_xy, strainrate_yy, ux_x, ux_y, uy_x, uy_y, mu, H)
    shearstress!(shearstress_x, shearstress_y, strainrate_xx, strainrate_xy, strainrate_yy, dx, dy)
    basalstress!(basalstress_x, basalstress_y, β_acx, β_acy, ux, uy)
    drivingstress!(drivingstress_x, drivingstress_y, rho_ice, g, H, z_b + H, dx, dy, nx, ny)
    pseudo_ux .= ( shearstress_x - basalstress_x ) ./ (rho_ice .* H) - drivingstress_x
    pseudo_uy .= ( shearstress_y - basalstress_y ) ./ (rho_ice .* H) - drivingstress_y
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
    (; mu) = state
    (; rho_ice, muB) = params
    (; dx, dy, nx, ny) = domain

    (; theta_v, theta_mu, reltol, maxiter) = options

    # Assign initial values
    pseudo_ux .= ux
    pseudo_uy .= uy
    pseudo_ux_old .= ux
    pseudo_uy_old .= uy

    # Init PT loop
    err = 2 * reltol
    iter = 0

    # PT loop following Sandip et al. (2024)
    while err > reltol || iter + 1 > maxiter
        pseudo_dotvel!(state, domain)
        dtau = pseudo_dt(rho, dx, mu, muB)
        pseudo_vel!(ux, ux_old, ux_dotvel, dtau, theta_v)
        pseudo_vel!(uy, uy_old, uy_dotvel, dtau, theta_v)
        err = norm(v - v_old)
        v_old .= v
        iter += 1
        @show err
    end
    
    return nothing
end
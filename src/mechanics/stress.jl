
"""
    calc_driving_stress!(icehseet)
    calc_driving_stress!(state, domain, params)
    calc_driving_stress!(taud_acx, taud_acy, H, z_srf, ρ, g, nx, ny, dx, dy)

Update the driving stress of the ice sheet.
"""
function calc_driving_stress!(icesheet::IceSheet)
    (; state, domain, params) = icesheet
    calc_driving_stress!(state, domain, params)
end

function calc_driving_stress!(state, domain, params)
    (; H, z_srf, taud_acx, taud_acy) = state
    (; nx, ny, dx, dy) = domain
    (; ρ, g) = params
    calc_driving_stress!(taud_acx, taud_acy, H, z_srf, ρ, g, nx, ny, dx, dy)
end

function calc_driving_stress!(taud_acx, taud_acy, H, z_srf, ρ, g, nx, ny, dx, dy)

    for i = 1:nx
        for j = 1:ny
            ip1 = periodic_bc_plusindices(i, nx)
            jp1 = periodic_bc_plusindices(j, ny)

            H_mid = 0.5 * (H[i, j] + H[ip1, j])
            taud_acx[i, j] = ρ * g * H_mid * (z_srf[ip1, j] - z_srf[i, j]) / dx

            H_mid = 0.5 * (H[i, j] + H[i, jp1])
            taud_acy[i, j] = ρ * g * H_mid * (z_srf[i, jp1] - z_srf[i, j]) / dy
        end
    end

    # BC: periodic...
    # Hack to ensure driving stress is ok in last grid point,
    # since z_srf might not be periodic...
    taud_acx[end, :] = taud_acx[end-1, :]
    return nothing
end


function shearstress!(shear_x, shear_y, strainrate_xx, strainrate_xy, strainrate_yy, dx, dy)
    shear_x .= delx(strainrate_xx, dx) + dely(strainrate_xy, dy)
    shear_y .= delx(strainrate_xy, dx) + dely(strainrate_yy, dy)
    return nothing
end

function basalstress!(basalstress_x, basalstress_y, β_acx, β_acy, u, v)
    basalstress_x .= β_acx .* u
    basalstress_y .= β_acy .* v
    return nothing
end

function drivingstress!(drivingstress_x, drivingstress_y, rho_ice, g, H, z_srf, dx, dy, nx, ny)
    drivingstress_x .= (rho_ice * g) .* H .* delx(z_srf, dx, nx)
    drivingstress_y .= (rho_ice * g) .* H .* dely(z_srf, dy, ny)
    return nothing
end
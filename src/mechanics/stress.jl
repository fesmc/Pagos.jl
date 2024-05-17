
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


function shearstress!(shear_x, shear_y, strainrate_xx, strainrate_xy, strainrate_yy,
    prealloc, dx, dy, nx, ny)
    
    # Allocation-free version of: shear_x .= delx(sr_xx, dx) + dely(sr_xy, dy)
    delx!(prealloc, strainrate_xx, dx, nx)
    shear_x .= prealloc
    dely!(prealloc, strainrate_xy, dy, ny)
    shear_x .+= prealloc

    # Allocation-free version of: shear_y .= delx(sr_xy, dx) + dely(sr_yy, dy)
    delx!(prealloc, strainrate_xy, dx, nx)
    shear_y .= prealloc
    dely!(prealloc, strainrate_yy, dy, ny)
    shear_y .+= prealloc
    return nothing
end

function basalstress!(basalstress_x, basalstress_y, β_acx, β_acy, u, v)
    basalstress_x .= β_acx .* u
    basalstress_y .= β_acy .* v
    return nothing
end

function drivingstress!(drivingstress_x, drivingstress_y, prealloc, rho_ice, g, H, z_b, dx, dy, nx, ny)
    @. prealloc = H + z_b
    delx!(drivingstress_x, prealloc, dx, nx)
    dely!(drivingstress_y, prealloc, dy, ny)
    drivingstress_x .*= (rho_ice * g) .* H
    drivingstress_y .*= (rho_ice * g) .* H
    return nothing
end
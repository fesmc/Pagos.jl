# TODO maybe need to stagger the stresses!
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

function basalstress!(basalstress_x, basalstress_y, β_acx, β_acy, ux, uy)
    basalstress_x .= β_acx .* ux
    basalstress_y .= β_acy .* uy
    return nothing
end

function drivingstress!(drivingstress_x, drivingstress_y, prealloc, rho_ice, g, H, z_b, dx, dy, nx, ny)
    @. prealloc = H + z_b
    delx!(drivingstress_x, prealloc, dx, nx)
    dely!(drivingstress_y, prealloc, dy, ny)
    @. drivingstress_x *= rho_ice * g * H
    @. drivingstress_y *= rho_ice * g * H
    return nothing
end
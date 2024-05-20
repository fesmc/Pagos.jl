"""
    shearstress!(shear_x, shear_y, strainrate_xx, strainrate_xy, strainrate_yy,
        prealloc, dx, dy, nx, ny)

Compute the shear stress components `shear_x` and `shear_y` from the components
of the scaled strain rate tensor, as computed by [`scaledstrainrate!`](@ref).
"""
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

"""
    basalstress!(basalstress_x, basalstress_y, beta_acx, beta_acy, ux, uy)

Compute the basal stress components `basalstress_x` and `basalstress_y` from the
basal friction coefficients `beta_acx` and `beta_acy` and the velocity components `ux` and `uy`.
"""
function basalstress!(basalstress_x, basalstress_y, beta_acx, beta_acy, ux, uy)
    basalstress_x .= beta_acx .* ux
    basalstress_y .= beta_acy .* uy
    return nothing
end

"""
    drivingstress!(drivingstress_x, drivingstress_y, prealloc, rho_ice, g, H, z_b, dx, dy, nx, ny)

Compute the driving stress components `drivingstress_x` and `drivingstress_y` from the
ice density `rho_ice`, the acceleration due to gravity `g`, the ice thickness `H`, the
bedrock elevation `z_b`, and the grid spacings `dx` and `dy`. The helper `prealloc` is
merely used for temporary storage.
"""
function drivingstress!(drivingstress_x, drivingstress_y, prealloc, rho_ice, g, H, z_b, dx, dy, nx, ny)
    @. prealloc = H + z_b
    delx!(drivingstress_x, prealloc, dx, nx)
    dely!(drivingstress_y, prealloc, dy, ny)
    @. drivingstress_x *= rho_ice * g * H
    @. drivingstress_y *= rho_ice * g * H
    return nothing
end
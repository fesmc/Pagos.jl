"""
$(TYPEDSIGNATURES)

Compute the deviatoric stress tensor of `mech` from its strain-rate tensor and the
ice viscosity `mat.eta_ice`, following the constitutive relation `` \\tau_{ij} = 2 \\, \\eta \\, \\dot\\varepsilon_{ij} ``,
and the effective (second-invariant) stress `mech.stress_effective`.

The computation is **independent of the dynamics**: the strain-rate, stress and
viscosity fields are all 3D `(nx, ny, nz)` arrays (depth-averaged solvers simply use
`nz == 1`), so every component is obtained the same way. In particular `` \\tau_{zz} ``
is reconstructed from the traceless property of the deviatoric stress tensor,
`` \\tau_{zz} = -(\\tau_{xx} + \\tau_{yy}) ``, which under incompressibility equals
`` 2\\eta\\dot\\varepsilon_{zz} `` exactly, so reconstructing loses no information. The
vertical-shear components `` \\tau_{xz}, \\tau_{yz} `` follow directly from
`` \\dot\\varepsilon_{xz}, \\dot\\varepsilon_{yz} ``, which are zero in the depth-averaged
case (no resolved vertical shear) and populated in the full-column case.

A **single fused KernelAbstractions kernel** does the whole tensor in one pass on
whichever backend (CPU, GPU) owns the `mech` arrays: it reads the viscosity and the
five independent strain-rate components once, then writes all six stress components
and the effective stress. Because every field is a 3D array of the same shape, the
kernel uses a flat `@index(Global, Linear)` index, so no 2D/3D special-casing is
needed; the single read/write sweep over memory is the bandwidth-optimal layout for
this otherwise memory-bound computation.

This assumes the `strain_rate_d*` fields of `mech` are already populated.
"""
function deviatoric_stress!(mech::MechanicState, mat::MaterialState)
    deviatoric_stress!(
        mech.stress.xx, mech.stress.yy, mech.stress.zz,
        mech.stress.xy, mech.stress.xz, mech.stress.yz, mech.stress.effective,
        mat.eta_ice,
        mech.strainrate.xx, mech.strainrate.yy,
        mech.strainrate.xy, mech.strainrate.xz, mech.strainrate.yz,
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

Array-level method of [`deviatoric_stress!`](@ref): launches the fused kernel directly
on the stress outputs (`sxx … seff`), the viscosity `η` and the five independent
strain-rate components (`exx, eyy, exy, exz, eyz`). All arrays must share the same shape
and reside on the same backend. `deviatoric_stress!(mech, mat)` is a thin wrapper that
unpacks the corresponding fields of `mech`/`mat` and calls this method.
"""
function deviatoric_stress!(sxx, syy, szz, sxy, sxz, syz, seff, η, exx, eyy, exy, exz, eyz)
    backend = get_backend(sxx)
    kernel! = _deviatoric_stress_kernel!(backend)
    kernel!(sxx, syy, szz, sxy, sxz, syz, seff, η, exx, eyy, exy, exz, eyz;
        ndrange = length(sxx))
    KernelAbstractions.synchronize(backend)
    return nothing
end
# TODO this typically does not need to be computed where we don't have any ice and could be easily handled via a mask passed to the kernel. Check performance!

@kernel function _deviatoric_stress_kernel!(
    sxx, syy, szz, sxy, sxz, syz, seff,
    η, exx, eyy, exy, exz, eyz,
)
    I = @index(Global, Linear)
    @inbounds begin
        twoη = 2 * η[I]
        txx  = twoη * exx[I]
        tyy  = twoη * eyy[I]
        txy  = twoη * exy[I]
        txz  = twoη * exz[I]
        tyz  = twoη * eyz[I]
        tzz  = -(txx + tyy)                      # traceless deviatoric identity
        sxx[I] = txx; syy[I] = tyy; szz[I] = tzz
        sxy[I] = txy; sxz[I] = txz; syz[I] = tyz
        seff[I] = sqrt((txx^2 + tyy^2 + tzz^2) / 2 + txy^2 + txz^2 + tyz^2)
    end
end

"""
    shearstress!(shear_x, shear_y, strainrate_xx, strainrate_xy, strainrate_yy,
        prealloc, dx, dy, nx, ny)

Compute the shear stress components `shear_x` and `shear_y` from the components
of the scaled strain rate tensor, as computed by [`scaledstrainrate!`](@ref).
"""
# TODO maybe need to stagger the stresses!
function shearstress!(shear_x, shear_y, strainrate_xx, strainrate_xy, strainrate_yy,
    prealloc, dx, dy)

    ∂x!(prealloc, strainrate_xx, dx)
    shear_x .= prealloc
    ∂y!(prealloc, strainrate_xy, dy)
    shear_x .+= prealloc

    ∂x!(prealloc, strainrate_xy, dx)
    shear_y .= prealloc
    ∂y!(prealloc, strainrate_yy, dy)
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
function drivingstress!(drivingstress_x, drivingstress_y, prealloc, rho_ice, g, H, z_b, dx, dy)
    @. prealloc = H + z_b
    ∂x!(drivingstress_x, prealloc, dx)
    ∂y!(drivingstress_y, prealloc, dy)
    @. drivingstress_x *= rho_ice * g * H
    @. drivingstress_y *= rho_ice * g * H
    return nothing
end
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

# TODO maybe need to stagger the stresses!
"""
    shearstress!(shear_x, shear_y, strainrate_xx, strainrate_xy, strainrate_yy,
        prealloc, dx, dy)

Compute the shear stress components `shear_x` and `shear_y` from the components
of the scaled strain rate tensor.
"""
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
$(TYPEDSIGNATURES)

Compute the basal-stress components `stress.base_x` and `stress.base_y` of `m` in place.

The basal stress follows the friction law `` \\boldsymbol{\\tau}_b = \\beta\\,\\mathbf{v}_b ``,
i.e. the effective basal friction coefficient `friction.beta_eff` times the basal velocity
`velocity.x_base`, `velocity.y_base`. Like the driving stress it is a depth-averaged 2D
field and dynamics-independent, so no dispatch on the momentum balance is needed.

Staggering of `beta_eff` onto the velocity points is assumed to be handled externally, so
the components are formed by a plain elementwise product.
"""
function basalstress!(m::Mechanics)
    (; state) = m
    (; stress, friction, velocity) = state
    return basalstress!(
        stress.base_x, stress.base_y,
        friction.beta_eff, velocity.x_base, velocity.y_base,
    )
end

"""
$(TYPEDSIGNATURES)

Array-level method of [`basalstress!`](@ref): writes the basal-stress components `base_x`,
`base_y` as the elementwise product of the effective basal friction coefficient `β` and the
basal velocity components `v_x`, `v_y`.
"""
function basalstress!(base_x, base_y, β, v_x, v_y)
    @. base_x = β * v_x
    @. base_y = β * v_y
    return nothing
end
# TODO: this should be staggered!

"""
$(TYPEDSIGNATURES)

Compute the driving-stress components `stress.driving_x` and `stress.driving_y` of `m`
in place.

The driving stress is `` \\tau_{d} = \\rho_{ice}\\,g\\,H\\,\\nabla s ``, where `` s `` is the
ice surface elevation, `` H `` the thickness and `` \\rho_{ice}, g `` are taken from the
physical `Constants`. It is **independent of the dynamics**: every momentum balance shares
this expression, so the depth-averaged 2D `driving_*` fields are obtained the same way and
no dispatch on the momentum balance is needed.
"""
function drivingstress!(m::Mechanics, c::Constants)
    (; state, grid) = m
    (; stress, topography) = state
    (; dx, dy) = grid
    return drivingstress!(
        stress.driving_x, stress.driving_y,
        topography.surface, topography.thickness,
        c.density_ice, c.gravity, dx, dy,
    )
end

"""
$(TYPEDSIGNATURES)

Array-level method of [`drivingstress!`](@ref): writes the driving-stress components
`driving_x`, `driving_y` from the surface elevation `surface`, the thickness `thickness`,
the ice density `ρ_ice`, the gravitational acceleration `g` and the grid spacings `dx`,
`dy`. The surface gradient is taken with the central-difference stencils `∂x!`/`∂y!`.
"""
function drivingstress!(driving_x, driving_y, surface, thickness, ρ_ice, g, dx, dy)
    ∂x!(driving_x, surface, dx)
    ∂y!(driving_y, surface, dy)
    @. driving_x *= ρ_ice * g * thickness
    @. driving_y *= ρ_ice * g * thickness
    return nothing
end
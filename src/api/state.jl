"""
    State{T<:AbstractFloat}

Struct containing the state variables of the ice sheet model.

# Fields
- `H`: ice thickness.
- `z_b`: bed elevation.
- `beta`: friction coefficient.
- `beta_acx`: friction coefficient staggered in x-direction.
- `beta_acy`: friction coefficient staggered in y-direction.
- `ux`: ice velocity in x.
- `uy`: ice velocity in y.
- `ux_old`: pseudo transient ice velocity in x at last iteration.
- `uy_old`: pseudo transient ice velocity in y at last iteration.
- `ux_x`: dux/dx.
- `ux_y`: dux/dy.
- `uy_x`: duy/dx.
- `uy_y`: duy/dy.
- `ux_b`: ux at bed.
- `uy_b`: uy at bed.
- `strainrate_xx`: xx-component of the strain rate tensor.
- `strainrate_xy`: xy-component of the strain rate tensor.
- `strainrate_yy`: yy-component of the strain rate tensor.
- `shearstress_x`: x-component of the shear stress.
- `shearstress_y`: y-component of the shear stress.
- `basalstress_x`: x-component of the basal stress.
- `basalstress_y`: y-component of the basal stress.
- `drivingstress_x`: x-component of the driving stress.
- `drivingstress_y`: y-component of the driving stress.
- `c_bed`: friction coefficient at bed.
- `f_ice`: fraction of ice.
- `mu`: ice viscosity.
- `N_ab`: effective viscosity.
- `prealloc`: temporary storage for calculations.
"""
mutable struct State{M}
    H::M
    z_b::M
    beta::M
    beta_acx::M
    beta_acy::M
    ux::M
    uy::M
    ux_old::M
    uy_old::M
    dotvel_x::M
    dotvel_y::M
    ux_x::M
    ux_y::M
    uy_x::M
    uy_y::M
    ux_b::M
    uy_b::M
    strainrate_xx::M
    strainrate_xy::M
    strainrate_yy::M
    shearstress_x::M
    shearstress_y::M
    basalstress_x::M
    basalstress_y::M
    drivingstress_x::M
    drivingstress_y::M
    c_bed::M
    f_ice::M
    mu::M
    N_ab::M
    prealloc::M
end
# TODO should use lazy maps instead of large structs

function State(domain::Domain{T}) where {T<:AbstractFloat}
    return State([copy(domain.null) for _ in fieldnames(State)]...)
end
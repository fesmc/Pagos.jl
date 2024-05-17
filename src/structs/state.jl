"""
    State{T<:AbstractFloat}

Struct containing the state variables of the ice sheet model, which contains:
- `H::Matrix{T}`: ice thickness.
- `z_b::Matrix{T}`: bed elevation.
- `taud_acx::Matrix{T}`: x-component of the shear stress.
- `taud_acy::Matrix{T}`: y-component of the shear stress.
- `β::Matrix{T}`: friction coefficient.
- `β_acx::Matrix{T}`: ?
- `β_acy::Matrix{T}`: ?
- `ux::Matrix{T}`: ice velocity in x.
- `uy::Matrix{T}`: ice velocity in y.
- `pseudo_ux::Matrix{T}`: pseudo transient ice velocity in x.
- `pseudo_uy::Matrix{T}`: pseudo transient ice velocity in y.
- `pseudo_ux_old::Matrix{T}`: pseudo transient ice velocity in x at last iteration.
- `pseudo_uy_old::Matrix{T}`: pseudo transient ice velocity in y at last iteration.
- `ux_x::Matrix{T}`: dux/dx.
- `ux_y::Matrix{T}`: dux/dy.
- `uy_x::Matrix{T}`: duy/dx.
- `uy_y::Matrix{T}`: duy/dy.
- `ux_b::Matrix{T}`: ux at bed.
- `uy_b::Matrix{T}`: uy at bed.
- `strainrate_xx::Matrix{T}`: xx-component of the strain rate tensor.
- `strainrate_xy::Matrix{T}`: xy-component of the strain rate tensor.
- `strainrate_yy::Matrix{T}`: yy-component of the strain rate tensor.
- `shearstress_x::Matrix{T}`: x-component of the shear stress.
- `shearstress_y::Matrix{T}`: y-component of the shear stress.
- `basalstress_x::Matrix{T}`: x-component of the basal stress.
- `basalstress_y::Matrix{T}`: y-component of the basal stress.
- `drivingstress_x::Matrix{T}`: x-component of the driving stress.
- `drivingstress_y::Matrix{T}`: y-component of the driving stress.
- `c_bed::Matrix{T}`: friction coefficient at bed.
- `f_ice::Matrix{T}`: ?
- `mu::Matrix{T}`: ice viscosity.
"""
mutable struct State{T<:AbstractFloat}
    H::Matrix{T}
    z_b::Matrix{T}
    taud_acx::Matrix{T}
    taud_acy::Matrix{T}
    β::Matrix{T}
    β_acx::Matrix{T}
    β_acy::Matrix{T}
    ux::Matrix{T}
    uy::Matrix{T}
    pseudo_ux::Matrix{T}
    pseudo_uy::Matrix{T}
    pseudo_ux_old::Matrix{T}
    pseudo_uy_old::Matrix{T}
    dotvel_x::Matrix{T}
    dotvel_y::Matrix{T}
    ux_x::Matrix{T}
    ux_y::Matrix{T}
    uy_x::Matrix{T}
    uy_y::Matrix{T}
    ux_b::Matrix{T}
    uy_b::Matrix{T}
    strainrate_xx::Matrix{T}
    strainrate_xy::Matrix{T}
    strainrate_yy::Matrix{T}
    shearstress_x::Matrix{T}
    shearstress_y::Matrix{T}
    basalstress_x::Matrix{T}
    basalstress_y::Matrix{T}
    drivingstress_x::Matrix{T}
    drivingstress_y::Matrix{T}
    c_bed::Matrix{T}
    f_ice::Matrix{T}
    mu::Matrix{T}
    prealloc::Matrix{T}
end
# TODO should use lazy maps instead of large structs

function State(domain::Domain{T}) where {T<:AbstractFloat}
    return State([copy(domain.null) for _ in fieldnames(State)]...)
end
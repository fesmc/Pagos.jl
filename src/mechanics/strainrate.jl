"""
$(TYPEDSIGNATURES)

Compute the strain rate tensor in place.
"""
function strainrate!(m::Mechanics)
    (; state, momentum) = m
    (; strainrate, velocity, material, topography) = state
    return strainrate!(
        strainrate, velocity, material, topography, momentum,
    )
end

function strainrate!(strainrate, velocity, material, topo, momentum::AbstractMomentumBalance)
    backend = get_backend(strainrate.xx)
    kernel! = _strainrate_kernel!(backend)
    kernel!(strainrate, velocity, material, topo, momentum; ndrange = length(strainrate.xx))
    KernelAbstractions.synchronize(backend)
    return nothing
end

# Per-element fallback: triggers when a momentum balance has no specialized method below,
# i.e. the real extension point. Guarding the launcher instead would never catch this.
function strainrate!(strainrate, velocity, material, topo, momentum::AbstractMomentumBalance, I)
    throw(ArgumentError("Unsupported momentum balance type: $(typeof(momentum))"))
end

function strainrate!(strainrate, velocity, material, topo, momentum::SIAMomentumBalance, I)
    # For SIA: depth-averaged viscosity, thickness, strainrate.xz/yz, velocity.*_bar_dz are all 2D
    strainrate.xz[I] = material.viscosity_depthaveraged[I] * topo.thickness[I] * velocity.x_bar_dz[I]
    strainrate.yz[I] = material.viscosity_depthaveraged[I] * topo.thickness[I] * velocity.y_bar_dz[I]
end

function strainrate!(strainrate, velocity, material, topo, momentum::MB, I) where MB<:Union{SSAMomentumBalance, DIVAMomentumBalance}
    # For SSA/DIVA: depth-averaged viscosity, thickness, strainrate.xx/xy/yy, velocity.*_d* are all 2D
    strainrate.xx[I] = 2 * material.viscosity_depthaveraged[I] * topo.thickness[I] * (2 * velocity.x_dx[I] + velocity.y_dy[I])
    strainrate.xy[I] = material.viscosity_depthaveraged[I] * topo.thickness[I] * (velocity.x_dy[I] + velocity.y_dx[I])
    strainrate.yy[I] = 2 * material.viscosity_depthaveraged[I] * topo.thickness[I] * (velocity.x_dx[I] + 2 * velocity.y_dy[I])
end

function strainrate!(strainrate, velocity, material, topo, momentum::BlatterPattynMomentumBalance, I)
    # For Blatter-Pattyn: 3D viscosity, strainrate.xx/xy/yy/xz/yz, velocity.*_d* are all 3D
    strainrate.xx[I] = 2 * material.viscosity[I] * (2 * velocity.x_dx[I] + velocity.y_dy[I])
    strainrate.xy[I] = material.viscosity[I] * (velocity.x_dy[I] + velocity.y_dx[I])
    strainrate.yy[I] = 2 * material.viscosity[I] * (velocity.x_dx[I] + 2 * velocity.y_dy[I])
    strainrate.xz[I] = material.viscosity[I] * velocity.x_dz[I]
    strainrate.yz[I] = material.viscosity[I] * velocity.y_dz[I]
end

@kernel function _strainrate_kernel!(strainrate, velocity, material, topo, momentum::AbstractMomentumBalance)
    I = @index(Global, Linear)
    @inbounds begin
        strainrate!(strainrate, velocity, material, topo, momentum, I)
        strainrate_effective!(strainrate, velocity, momentum, I) 
    end
end

"""
$(TYPEDSIGNATURES)

Compute the effective strain rate in place.
"""
# Per-element fallback: triggers when a momentum balance has no specialized method below.
function strainrate_effective!(strainrate, velocity, momentum::AbstractMomentumBalance, I)
    throw(ArgumentError("Unsupported momentum balance type: $(typeof(momentum))"))
end

function strainrate_effective!(strainrate, velocity, momentum::SIAMomentumBalance, I)
    strainrate.effective[I] = sqrt(1 / 4 * (velocity.x_bar_dz[I] + velocity.y_bar_dz[I]) ^ 2)
end

function strainrate_effective!(strainrate, velocity, momentum::SSAMomentumBalance, I)
    strainrate.effective[I] = sqrt(
        velocity.x_dx[I]^2 + velocity.y_dy[I]^2 +
        velocity.x_dx[I] * velocity.y_dy[I] +
        1 / 4 * (velocity.x_dy[I] + velocity.y_dx[I]) ^ 2
    )
end

function strainrate_effective!(strainrate, velocity, momentum::MB, I) where MB<:Union{DIVAMomentumBalance, BlatterPattynMomentumBalance}
    strainrate.effective[I] = sqrt(
        velocity.x_dx[I]^2 + velocity.y_dy[I]^2 +
        velocity.x_dx[I] * velocity.y_dy[I] +
        1 / 4 * (velocity.x_dy[I] + velocity.y_dx[I]) ^ 2 +
        1 / 4 * velocity.x_dz[I]^2 +
        1 / 4 * velocity.y_dz[I]^2
    )
end


"""
$(TYPEDSIGNATURES)

Compute the velocity gradients in x (`v_x_dx, v_y_x`) and y-direction (`v_x_dy, v_y_y`).
The gradients are computed using the central difference scheme. The input velocities
`v_x` and `v_y` are defined on a staggered grid with dimensions `nx` and `ny`.
The grid spacing in x and y-direction is given by `dx` and `dy`.
"""
function velocitygradients!(mechanics::Mechanics)
    (; state, grid) = mechanics
    (; velocity) = state
    (; dx, dy) = grid
    return velocitygradients!(velocity, dx, dy)
end

function velocitygradients!(velocity::VelocityState, dx, dy)
    (; x, y) = velocity
    return velocitygradients!(velocity.x_dx, velocity.x_dy, velocity.y_dx, velocity.y_dy, x, y, dx, dy)
end

function velocitygradients!(v_x_dx, v_x_dy, v_y_dx, v_y_dy, v_dx, v_dy, dx, dy)
    ∂x!(v_x_dx, v_dx, dx)
    ∂y!(v_x_dy, v_dx, dy)
    ∂x!(v_y_dx, v_dy, dx)
    ∂y!(v_y_dy, v_dy, dy)
    return nothing
end
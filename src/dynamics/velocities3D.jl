"""
    aggregate_viscosity_integral!(Fm, mu, H, m, sigma, l)

Aggregate the viscosity integral for the `l`-th layer to the final result `Fm`.

# Arguments
- `Fm::Matrix{T}`: Final viscosity integral.
- `mu::Array{T, 3}`: Viscosity tensor.
- `H::Matrix{T}`: Ice thickness.
- `m::Int`: Exponent.
- `sigma::Vector{T}`: Vertical (transformed) coordinate.
- `l::Int`: Layer index.
"""
# TODO: we should apply a nicer numerical integration scheme here.
# The current implementation is a simple Riemann sum. Gauss-Legendre quadrature
# is already available but would imply the assumption of a linear interpolation
# over z. This might slow down the computation.
function aggregate_viscosity_integral!(
    Fm::Matrix{T},
    mu::Array{T, 3},
    H::Matrix{T},
    m::Int,
    sigma::Vector{T},
    l::Int,
) where {T<:AbstractFloat}
    s_minus_z_over_H = 1-sigma[l]
    if l == 1
        dsigma = sigma[l]
    else
        dsigma = sigma[l] - sigma[l-1]
    end
    @. Fm += ( s_minus_z_over_H ^ m * view(dsigma, l) * H ) / view(mu, :, :, l)
end

"""
    aggregated_viscosity_integral!(Fm, mu, H, m, sigma)

Aggregate the viscosity integral for all layers to the final result `Fm`.

# Arguments
- `Fm::Matrix{T}`: Final viscosity integral.
- `mu::Array{T, 3}`: Viscosity tensor.
- `H::Matrix{T}`: Ice thickness.
- `m::Int`: Exponent.
- `sigma::Vector{T}`: Vertical (transformed) coordinate.
"""
function aggregated_viscosity_integral!(
    Fm::Matrix{T},
    mu::Array{T, 3},
    H::Matrix{T},
    m::Int,
    sigma::Vector{T},
) where {T<:AbstractFloat}

    Fm .= T(0)
    for l in eachindex(sigma)
        aggregate_viscosity_integral!(Fm, mu, H, m, sigma, l)
    end
    return nothing
end

"""
    layer_velocity(ub, beta, H, F1)

Compute the 2D velocity field for a given layer, which is represented by the `F1` integral.

# Arguments
- `ub::Matrix{T}`: Basal velocity.
- `beta::Matrix{T}`: Friction coefficient.
- `H::Matrix{T}`: Ice thickness.
- `F1::Matrix{T}`: Generalized viscosity integral at layer `l`.
"""
layer_velocity(ub, beta, H, F1) = ub .+ beta .* ub .* F1 ./ H

"""
    velocities3D!(F1, vx, vy, vz, mu, beta, H, sigma)

Compute the 3D velocity field.

# Arguments
- `F1::Matrix{T}`: Generalized viscosity integral.
- `vx::Array{T, 3}`: x-component of the velocity field.
- `vy::Array{T, 3}`: y-component of the velocity field.
- `vz::Array{T, 3}`: z-component of the velocity field.
- `mu::Array{T, 3}`: Viscosity tensor.
- `beta::Matrix{T}`: Friction coefficient.
- `H::Matrix{T}`: Ice thickness.
- `sigma::Vector{T}`: Vertical (transformed) coordinate.
"""
function velocities3D!(
    F1::Matrix{T},
    vx::Array{T, 3},
    vy::Array{T, 3},
    vz::Array{T, 3},
    mu::Array{T, 3},
    beta::Matrix{T},
    H::Matrix{T},
    sigma::Vector{T},
) where {T<:AbstractFloat}

    F1 .= T(0)
    for l in eachindex(sigma)
        aggregate_viscosity_integral!(F1, mu, H, 1, sigma, l)
        vx[:, :, l] .= layer_velocity(vx[:, :, l], beta, H, F1)
        vy[:, :, l] .= layer_velocity(vy[:, :, l], beta, H, F1)
        # vz[:, :, l] .= tau_b_x * (s - z) / (eta(z) * H)
    end
    
    return nothing
end


"""
    surface_velocity!(u_s, u_b, beta, F1)
    surface_velocity(u_b, beta, F1)

Compute the surface velocity `u_s` from the basal velocity `u_b` and the friction
coefficient `beta`. If possible, use the in-place version `surface_velocity!` to
avoid unnecessary memory allocations.

# Arguments
- `u_s::Matrix{T}`: Surface velocity.
- `u_b::Matrix{T}`: Basal velocity.
- `beta::Matrix{T}`: Friction coefficient.
- `F1::Matrix{T}`: Generalized viscosity integral.
"""

function surface_velocity!(
    u_s::Matrix{T},
    u_b::Matrix{T},
    beta::Matrix{T},
    F1::Matrix{T},
) where {T<:AbstractFloat}
    @. u_s = u_b * ( 1 + beta * F1 )
    return nothing
end

function surface_velocity(
    u_b::Matrix{T},
    beta::Matrix{T},
    F1::Matrix{T},
) where {T<:AbstractFloat}

    nx, ny = size(u_b)
    u_s = zeros(T, nx, ny)
    surface_velocity!(u_s, u_b, beta, F1)

    return u_s
end

"""

    depthaveraged_velocity!(u_d, u_b, beta, F2)
    depthaveraged_velocity(u_b, beta, F2)

Compute the depth-averaged velocity `u_d` from the basal velocity `u_b` and the friction
coefficient `beta`. If possible, use the in-place version `depthaveraged_velocity!` to
avoid unnecessary memory allocations.

# Arguments
- `u_d::Matrix{T}`: Depth-averaged velocity.
- `u_b::Matrix{T}`: Basal velocity.
- `beta::Matrix{T}`: Friction coefficient.
- `F2::Matrix{T}`: Generalized viscosity integral.
"""
function depthaveraged_velocity!(
    u_d::Matrix{T},
    u_b::Matrix{T},
    beta::Matrix{T},
    F2::Matrix{T},
) where {T<:AbstractFloat}
    @. u_d = u_b * ( 1 + beta * F2 )
    return nothing
end

function depthaveraged_velocity(
    u_b::Matrix{T},
    beta::Matrix{T},
    F2::Matrix{T},
) where {T<:AbstractFloat}

    nx, ny = size(u_b)
    u_d = zeros(T, nx, ny)
    depthaveraged_velocity!(u_d, u_b, beta, F2)

    return u_d
end
























# """
#     generalized_viscosity_integral(mu, H, s, m, z, dz)
#     generalized_viscosity_integral!(Fm, mu, H, s, m, z, dz)

# Compute the generalized viscosity integral `Fm`, as introduced in Arthern et al. (2015).
# If possible, use the in-place version `generalized_viscosity_integral!` to avoid unnecessary
# memory allocations.

# # Arguments
# - `mu::Array{T, 3}`: Viscosity tensor.
# - `H::Matrix{T}`: Ice thickness.
# - `s::Matrix{T}`: Surface elevation.
# - `m::Int`: Exponent.
# - `z::Vector{T}`: Vertical coordinate.
# - `dz::Vector{T}`: Vertical grid spacing.

# Note: make sure that the vertical coordinate is not transformed!
# """
# function generalized_viscosity_integral!(
#     Fm::Matrix{T},
#     mu::Array{T, 3},
#     H::Matrix{T},
#     s::Matrix{T},
#     m::Int,
#     z::Vector{T},
#     dz::Vector{T},
# ) where {T<:AbstractFloat}

#     Fm .= T(0)
#     for l in eachindex(dz)
#         @. Fm += ((s - view(z, l)) / H) ^ m / mu[:, :, l] * dz[l]
#     end
#     return nothing
# end


# function generalized_viscosity_integral(
#     mu::Array{T, 3},
#     H::Matrix{T},
#     s::Matrix{T},
#     m::Int,
#     z::Vector{T},
#     dz::Vector{T},
# ) where {T<:AbstractFloat}

#     nx, ny = size(H)
#     Fm = zeros(T, nx, ny)
#     generalized_viscosity_integral!(Fm, mu, H, s, m, z, dz)

#     return Fm
# end





# # Functions to calculate velocity 
# function F_integral!(F, mu, H, f, zeta_aa, m, nx, ny)
#     for i = 1:nx, j = 1:ny
#         F[i, j] = 1 / mu[i, j] * H[i, j]^f * zeta_aa[i, j, k]
#     end
#     return Fn
# end

# function basalvelocity!(state, domain, params, tools)
#     (; H, z_b, ux, uy, beta_acx, beta_acy) = state
#     (; nx, ny, dx, dy) = domain
#     (; ρ, g) = params
#     (; Ai, Aj, Av, u, b) = tools
#     basalvelocity!(ux, uy, H, z_b, beta_acx, beta_acy, ρ, g, nx, ny, dx, dy, Ai, Aj, Av, u, b)
#     return nothing
# end
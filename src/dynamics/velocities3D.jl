"""
    basal_velocity_from_depthavg_velocity!(vb, v, beta, F2, mask_incides)

Compute the basal velocity `vb` from the depth-avergaed velocity field `v`.

# Arguments
- `vb::Matrix{T}`: Basal velocity.
- `v::Matrix{T}`: Velocity field.
- `beta::Matrix{T}`: Friction coefficient.
- `F2::Matrix{T}`: Generalized viscosity integral with `m=2`.
- `mask_incides::Vector{CartesianIndex{2}}`: Mask indices.
"""
function basal_velocity_from_depthavg_velocity!(
    vb::Matrix{T},
    v::Matrix{T},
    beta::Matrix{T},
    F2::Matrix{T},
    mask_incides::Vector{CartesianIndex{2}},
) where {T<:AbstractFloat}
    for I in mask_incides
        vb[I] = v[I] / (1 + beta[I] * F2[I])
    end
    return nothing
end

"""
    basal_velocity_from_surface_velocity!(vb, v, beta, F2, mask_incides)

Compute the basal velocity `vb` from the surface velocity `v`.

# Arguments
- `vb::Matrix{T}`: Basal velocity.
- `v::Matrix{T}`: Surface velocity.
- `beta::Matrix{T}`: Friction coefficient.
- `F2::Matrix{T}`: Generalized viscosity integral with `m=2`.
- `mask_incides::Vector{CartesianIndex{2}}`: Mask indices.
"""
function basal_velocity_from_surface_velocity!(
    vb::Matrix{T},
    vs::Matrix{T},
    beta::Matrix{T},
    F1::Matrix{T},
    mask_incides::Vector{CartesianIndex{2}},
) where {T<:AbstractFloat}
    for I in mask_incides
        vb[I] = vs[I] / (1 + beta[I] * F1[I])
    end
    return nothing
end

"""
    surface_velocity!(vs, vb, beta, F1, mask_incides)

Compute the surface velocity `vs` from the basal velocity `vb`.

# Arguments
- `vs::Matrix{T}`: Surface velocity.
- `vb::Matrix{T}`: Basal velocity.
- `beta::Matrix{T}`: Friction coefficient.
- `F1::Matrix{T}`: Generalized viscosity integral with `m=1`.
- `mask_incides::Vector{CartesianIndex{2}}`: Mask indices.
"""
function surface_velocity!(
    vs::Matrix{T},
    vb::Matrix{T},
    beta::Matrix{T},
    F1::Matrix{T},
    mask_incides::Vector{CartesianIndex{2}},
) where {T<:AbstractFloat}
    for I in mask_incides
        vs[I] = vb[I] * (1 + beta[I] * F1[I])
    end
    return nothing
end

"""
    depthavg_velocity!(v, vb, beta, F2, mask_incides)

Compute the depth-averaged velocity `v` from the basal velocity `vb`.

# Arguments
- `v::Matrix{T}`: Depth-averaged velocity.
- `vb::Matrix{T}`: Basal velocity.
- `beta::Matrix{T}`: Friction coefficient.
- `F2::Matrix{T}`: Generalized viscosity integral with `m=2`.
- `mask_incides::Vector{CartesianIndex{2}}`: Mask indices.
"""
function depthavg_velocity!(
    v::Matrix{T},
    vb::Matrix{T},
    beta::Matrix{T},
    F2::Matrix{T},
    mask_incides::Vector{CartesianIndex{2}},
) where {T<:AbstractFloat}
    for I in mask_incides
        v[I] = vb[I] * (1 + beta[I] * F2[I])
    end
    return nothing
end

"""
    velocities3D!(F1, vx3D, vy, vz, mu, beta, H, sigma)

Compute the 3D velocity field.

# Arguments
- `F1::Matrix{T}`: Generalized viscosity integral with `m=1`.
- `vx3D::Array{T, 3}`: x-component of the velocity field.
- `vy::Array{T, 3}`: y-component of the velocity field.
- `vz::Array{T, 3}`: z-component of the velocity field.
- `mu::Array{T, 3}`: Viscosity tensor.
- `beta::Matrix{T}`: Friction coefficient.
- `H::Matrix{T}`: Ice thickness.
- `sigma::Vector{T}`: Vertical (transformed) coordinate.
"""
function velocities3D!(
    F1::Matrix{T},
    vx3D::Array{T, 3},
    vy3D::Array{T, 3},
    vz3D::Array{T, 3},
    vb_x::Matrix{T},
    vb_y::Matrix{T},
    mu::Array{T, 3},
    beta::Matrix{T},
    H::Matrix{T},
    sigma::Vector{T},
    mask_incides::Vector{CartesianIndex{2}},
) where {T<:AbstractFloat}

    F1 .= T(0)
    for l in eachindex(sigma)
        aggregate_viscosity_integral!(F1, mu, H, 1, sigma, l, mask_incides)
        layer_velocity!(vx3D, l, vb_x, beta, H, F1, mask_incides)
        layer_velocity!(vy3D, l, vb_y, beta, H, F1, mask_incides)
        # vz[:, :, l] .= tau_b_x * (s - z) / (eta(z) * H)
    end
    
    return nothing
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
- `mask_incides::Vector{CartesianIndex{2}}`: Mask indices.
"""
function aggregated_viscosity_integral!(
    Fm::Matrix{T},
    mu::Array{T, 3},
    H::Matrix{T},
    m::Int,
    sigma::Vector{T},
    mask_incides::Vector{CartesianIndex{2}},
) where {T<:AbstractFloat}

    Fm .= T(0)
    for l in eachindex(sigma)
        aggregate_viscosity_integral!(Fm, mu, H, m, sigma, l, mask_incides)
    end
    return nothing
end

"""
    aggregate_viscosity_integral!(Fm, mu, H, m, sigma, l, mask_indices)

Aggregate the viscosity integral for the `l`-th layer to the final result `Fm`.

# Arguments
- `Fm::Matrix{T}`: Final viscosity integral.
- `mu::Array{T, 3}`: Viscosity tensor.
- `H::Matrix{T}`: Ice thickness.
- `m::Int`: Exponent.
- `sigma::Vector{T}`: Vertical (transformed) coordinate.
- `l::Int`: Layer index.
- `mask_incides::Vector{CartesianIndex{2}}`: Mask indices.
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
    mask_incides::Vector{CartesianIndex{2}},
) where {T<:AbstractFloat}
    if l == 1
        dsigma = sigma[l]
    else
        dsigma = sigma[l] - sigma[l-1]
    end

    s_minus_z_over_H = 1-sigma[l]+dsigma/2

    for I in mask_incides
        Fm[I] += ( s_minus_z_over_H ^ m * dsigma * H[I] ) / mu[I, l]    # dsigma * H = dz
    end
end

"""
    layer_velocity!(v, l, vb, beta, H, F1, mask_incides)

Compute the velocity field for a given layer, which is represented by the `F1` integral.
To obtain the velocity field in `x` and `y` direction, provide the basal velocity `vb_x`
and `vb_y` respectively.

# Arguments
- `v::Array{T, 3}`: 3D velocity field.
- `l::Int`: Layer index.
- `vb::Matrix{T}`: Basal velocity.
- `beta::Matrix{T}`: Friction coefficient.
- `H::Matrix{T}`: Ice thickness.
- `F1::Matrix{T}`: Generalized viscosity integral at layer `l` with `m=1`.
- `mask_incides::Vector{CartesianIndex{2}}`: Mask indices.
"""
function layer_velocity!(
    v::Array{T, 3},
    l::Int,
    vb::Matrix{T},
    beta::Matrix{T},
    H::Matrix{T},
    F1::Matrix{T},
    mask_incides::Vector{CartesianIndex{2}},
) where {T<:AbstractFloat}
    for I in mask_incides
        v[I, l] = vb[I] + beta[I] * vb[I] * F1[I]
    end
    return nothing
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
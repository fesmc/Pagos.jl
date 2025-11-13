"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the basal friction computation via [`basal_shear_stress`](@ref).
"""
abstract type AbstractBasalFriction{T<:AbstractFloat} end

"""
$(TYPEDSIGNATURES)

A struct that defines a constant basal friction coefficient (beta). When passed to [`basal_shear_stress`](@ref), the basal shear stress is calculated as:

```math
\\begin{aligned}
\\tau_b = - \\beta \\, v_b
\\end{aligned}
```

# Fields:
 - `beta::T`: The constant basal friction coefficient.
"""
@kwdef struct ConstantBetaBasalFriction{T} <: AbstractBasalFriction{T}
    beta::T
end

"""
$(TYPEDSIGNATURES)

A struct that defines a linear relationship between basal friction coefficient (beta) and basal velocity.
"""
@kwdef struct LinearBetaBasalFriction{T} <: AbstractBasalFriction{T}
    u0::T
end

"""
$(TYPEDSIGNATURES)

A struct that defines the pseudo-plastic power law for basal sliding.
When passed to [`basal_shear_stress!`](@ref), the basal shear stress is calculated
following Eq. (24) of [robinson_description_2020](@citet). This covers plastic friction
(q = 0), linear friction (q = 1), and power-law friction (q > 1).
"""
@kwdef struct PseudoPlasticPowerBasalFriction{T} <: AbstractBasalFriction{T}
    v_0::T = 1e2        # "m/yr" reference velocity
    q::T = 0.8          # exponent
end

"""
    RegularizedCoulombBasalFriction{T}

A struct that defines the regularized Coulomb friction law for basal sliding.
When passed to [`basal_shear_stress!`](@ref), the basal shear stress is calculated
following Eq. (25) of [robinson_description_2020](@citet).
"""
@kwdef struct RegularizedCoulombBasalFriction{T} <: AbstractBasalFriction{T}
    v_0::T = 1e2        # "m/yr" reference velocity
    q::T = 0.2          # exponent
end

#=
"""
    AnisotropicRegularizedCoulombFriction{T}

Similar to [`RegularizedCoulombBasalFriction`](@ref), but with anisotropic friction coefficients
that can be arising through anisotropic bed roughness below resolution.
"""
struct AnisotropicRegularizedCoulombBasalFriction{T} <: AbstractBasalFriction{T}
    v_0::T
    q::T
end
=#

"""
$(TYPEDSIGNATURES)

Compute the basal shear stress `τ_basal` based on the basal velocity `v_basal`, the basal friction coefficient `c_basal`, and the basal friction model `friction<:AbstractBasalFriction`.
"""
function basal_shear_stress(
    v_basal::V,
    c_basal,
    friction::PseudoPlasticPowerBasalFriction,
) where {V<:AbstractVector{<:Real}}
    (; v_0, q) = friction
    v_basal_norm = norm(v_basal)
    return - c_basal * (v_basal_norm / v_0) ^ q * v_basal / v_basal_norm
end

function basal_shear_stress(
    v_basal::V,
    c_basal,
    friction::RegularizedCoulombBasalFriction,
) where {V<:AbstractVector{<:Real}}
    (; v_0, q) = friction
    v_basal_norm = norm(v_basal)
    return - c_basal * (v_basal_norm / (v_basal_norm + v_0)) ^ q * v_basal / v_basal_norm
end

function basal_shear_stress(
    v_basal,
    c_basal,
    friction,
)
    τ_basal = similar(v_basal)
    basal_shear_stress!(τ_basal, v_basal, c_basal, friction)
    return τ_basal
end

"""
$(TYPEDSIGNATURES)

Same as [`basal_shear_stress`](@ref) but operates in place.
"""
function basal_shear_stress!(
    τ_basal,
    v_basal,
    c_basal,
    friction::AbstractBasalFriction,
)
    map!((v, c) -> basal_shear_stress(v, c, friction), τ_basal, v_basal, c_basal)
    return nothing
end

"""
$(TYPEDSIGNATURES)

This will implement Eq. 2 of Zoet & Iverson (2020).
"""
function ut_with_clast()
    error("Not implemented yet")
end
#=
function basal_shear_stress!(
    v_basal_norm::Matrix{T},
    c_basal::Matrix{T},
    stress_basal_x::Matrix{T},
    stress_basal_y::Matrix{T},
    v_basal_x::Matrix{T},
    v_basal_y::Matrix{T},
    friction::RegularizedCoulombBasalFriction{T},
    mask::AbstractBitMask,
) where {T<:AbstractFloat}
    coulomb_basal_shear!(v_basal_norm, stress_basal_x, stress_basal_y, v_basal_x, v_basal_y,
        friction, c_basal, c_basal, mask)
    return nothing
end

function basal_shear_stress!(
    v_basal_norm::Matrix{T},
    c_basal_x::Matrix{T},
    c_basal_y::Matrix{T},
    stress_basal_x::Matrix{T},
    stress_basal_y::Matrix{T},
    v_basal_x::Matrix{T},
    v_basal_y::Matrix{T},
    friction::AnisotropicRegularizedCoulombFriction{T},
    mask::AbstractBitMask,
) where {T<:AbstractFloat}
    coulomb_basal_shear!(v_basal_norm, stress_basal_x, stress_basal_y, v_basal_x, v_basal_y,
        friction, c_basal_x, c_basal_y, mask)
    return nothing
end

function coulomb_basal_shear!(
    v_basal_norm::Matrix{T},
    stress_basal_x::Matrix{T},
    stress_basal_y::Matrix{T},
    v_basal_x::Matrix{T},
    v_basal_y::Matrix{T},
    friction::AbstractBasalFriction{T},
    c_basal_x::Matrix{T},
    c_basal_y::Matrix{T},
    mask::AbstractBitMask,
) where {T<:AbstractFloat}
    (; v_0, q) = friction
    @inbounds for I in view(mask)
        v_basal_norm[I] = norm(v_basal_x[I], v_basal_y[I])
        stress_basal_x[I] = coulomb_basal_shear(c_basal_x[I], v_basal_x[I],
            v_basal_norm[I], v_0, q)
        stress_basal_y[I] = coulomb_basal_shear(c_basal_y[I], v_basal_y[I],
            v_basal_norm[I], v_0, q)
    end
end

coulomb_basal_shear(c_basal, v_basal, v_basal_norm, v_0, q) = 
    - c_basal * (v_basal_norm / (v_basal_norm + v_0)) ^ q * v_basal / v_basal_norm

###############################################


"""
    calc_beta_aa_power_plastic(ux_b,uy_b,c_bed,f_ice,q,u_0)

Calculate basal friction coefficient (beta) that enters the SSA solver as a function
of basal velocity using a power-law form following Bueler and van Pelt (2015).    
"""
function calc_beta_aa_power_plastic(
    ux_b::Matrix{T},
    uy_b::Matrix{T},
    c_bed::Matrix{T},
    f_ice::Matrix{T},
    q,
    u_0,
) where {T}

    # Local variables
    ub_min = 1e-3               # [m/yr] Small min. velocity > 0 to avoid divide by 0
    ub_sq_min = ub_min^2
    nx, ny = size(ux_b)
    beta = fill(0.0, nx, ny)    # Initially set friction to zero everywhere

    for i = 1:nx
        for j = 1:ny
            im1, jm1 = periodic_minusindex(i, nx), periodic_minusindex(j, ny)

            if f_ice[i, j] == 1.0
                # Fully ice-covered point with some fully ice-covered neighbors 
                cb_aa = c_bed[i, j]

                if q == 1.0
                    # Linear law, no f(ub) term
                    beta[i, j] = cb_aa / u_0
                else
                    # Non-linear law with f(ub) term 
                    # Unstagger velocity components to aa-nodes 
                    ux_aa = 0.5 * (ux_b[i, j] + ux_b[im1, j])
                    uy_aa = 0.5 * (uy_b[i, j] + uy_b[i, jm1])
                    uxy_aa = sqrt(ux_aa^2 + uy_aa^2 + ub_sq_min)

                    if q == 0
                        # Plastic law
                        beta[i, j] = cb_aa * (1.0 / uxy_aa)
                    else
                        beta[i, j] = cb_aa * (uxy_aa / u_0)^q * (1.0 / uxy_aa)
                    end
                end

            else
                # Assign minimum velocity value, no staggering for simplicity

                if q == 1.0
                    # Linear law, no f(ub) term
                    beta[i, j] = c_bed[i, j] / u_0

                else
                    uxy_b = ub_min

                    if q == 0.0
                        # Plastic law
                        beta[i, j] = c_bed[i, j] * (1.0 / uxy_b)
                    else
                        beta[i, j] = c_bed[i, j] * (uxy_b / u_0)^q * (1.0 / uxy_b)
                    end
                end
            end
        end
    end

    return beta
end

function calc_beta_aa_power_plastic_nodes(
    ux_b::Matrix{T},
    uy_b::Matrix{T},
    c_bed::Matrix{T},
    f_ice::Matrix{T},
    q,
    u_0,
) where {T}


    # Local variables
    ub_min = 1e-3               # [m/yr] Minimum velocity is positive small value to avoid divide by zero
    ub_sq_min = ub_min^2

    nx, ny = size(ux_b)

    # Initially set friction to zero everywhere
    beta = fill(0.0, nx, ny)

    wt0 = 1.0 / sqrt(3)
    xn = [wt0, -wt0, -wt0, wt0]
    yn = [wt0, wt0, -wt0, -wt0]
    wtn = [1.0, 1.0, 1.0, 1.0]

    for i = 1:nx
        for j = 1:ny

            if f_ice[i, j] == 1.0
                # Fully ice-covered point with some fully ice-covered neighbors 
                cb_aa = c_bed[i, j]

                if q == 1.0
                    # Linear law, no f(ub) term
                    beta[i, j] = cb_aa / u_0

                else
                    # Non-linear law with f(ub) term 
                    uxn = acx_to_nodes(ux_b, i, j, xn, yn)
                    uyn = acy_to_nodes(uy_b, i, j, xn, yn)
                    uxyn = sqrt.(uxn .^ 2 .+ uyn .^ 2 .+ ub_sq_min)

                    if q == 0
                        # Plastic law
                        betan = cb_aa .* (1.0 ./ uxyn)
                    else
                        betan = cb_aa .* (uxyn ./ u_0) .^ q .* (1.0 ./ uxyn)
                    end
                    beta[i, j] = sum(betan .* wtn) / sum(wtn)
                end

            else
                # Assign minimum velocity value, no staggering for simplicity
                if q == 1.0
                    # Linear law, no f(ub) term
                    beta[i, j] = c_bed[i, j] / u_0
                else
                    uxy_b = ub_min
                    if q == 0.0
                        # Plastic law
                        beta[i, j] = c_bed[i, j] * (1.0 / uxy_b)
                    else
                        beta[i, j] = c_bed[i, j] * (uxy_b / u_0)^q * (1.0 / uxy_b)
                    end
                end
            end
        end
    end

    return beta
end
=#
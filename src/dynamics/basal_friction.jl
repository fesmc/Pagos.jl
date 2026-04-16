"""
$(TYPEDSIGNATURES)

A struct that encapsulates all basal friction related models and allows for easy update of the basal friction via [`basal_friction!`](@ref).

# Fields:
 - `beta<:AbstractBasalBeta`: The model for the basal friction coefficient.
 - `beta_gz<:AbstractBasalBetaGroundingZone`: The model for the basal friction coefficient at the grounding zone.
 - `roughness_sampling<:AbstractBedRoughnessSampling`: The model for the bed roughness sampling.
 - `c_bed_ref<:AbstractCbedRef`: The model for the reference bed friction coefficient.
 - `c_bed<:AbstractCbed`: The model for the bed friction coefficient.

```math
\\begin{aligned}
\\tau_b = - \\beta \\, v_b
\\end{aligned}
```
"""
struct BasalFriction{
    BB,     # <: AbstractBasalBeta
    BBGZ,   # <: AbstractBasalBetaGroundingZone
    BRS,    # <: AbstractBedRoughnessSampling
    CBR,    # <: AbstractCbedRef
    CB,     # <: AbstractCbed
}
    beta::BB
    beta_gz::BBGZ
    roughness_sampling::BRS
    c_bed_ref::CBR
    c_bed::CB
end

"""
$(TYPEDSIGNATURES)
"""
function basal_friction!(dyn_now, dyn_ref, topo_now, bf::BasalFriction)
    (; τ_basal, v_basal, β_basal, c_bed, N_eff) = dyn_now
    (; c_bed_ref) = dyn_ref
    (; z_bed, z_bed_σ, z_sl) = topo_now

end

function basal_friction!(
    τ_basal,
    v_basal,
    β_basal,
    c_bed,
    c_bed_ref,
    N_eff,
    z_bed,
    z_bed_σ,
    z_sl,
    bf::BasalFriction,
)

    c_bed!(c_bed, c_bed_ref, N_eff, bf.c_bed)
    basal_beta!(β_basal, v_basal, c_bed, c_bed_ref, bf.beta)
    basal_beta_gz!(β_basal, bf.beta_gz)
    basal_shear_stress!(τ_basal, v_basal, β_basal)
    return nothing
end


"""
$(TYPEDSIGNATURES)
"""
function basal_shear_stress(v_basal, β_basal)
    return - β_basal * v_basal
end
function basal_shear_stress(v_basal::M, β_basal::M) where {M<:AbstractMatrix{<:Real}}
    τ_basal = similar(v_basal)
    basal_shear_stress!(τ_basal, v_basal, β_basal)
    return τ_basal
end

"""
$(TYPEDSIGNATURES)
"""
function basal_shear_stress!(τ_basal, v_basal, β_basal)
    map!((v, β) -> basal_shear_stress(v, β), τ_basal, v_basal, β_basal)
    return nothing
end

#####################################################################
# Basal friction coefficient β
#####################################################################
"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the basal friction coefficient computation via [`basal_beta`](@ref).
"""
abstract type AbstractBasalBeta end

"""
$(TYPEDSIGNATURES)

A struct that defines a constant basal friction coefficient ``\beta``. When passed to [`basal_beta`](@ref), the basal shear stress is calculated as:

```math
\\begin{aligned}
\\beta = \\mathrm{const.}
\\end{aligned}
```

# Fields:
 - `beta::M`: The constant basal friction coefficient. Can be a scalar or an array.
"""
struct ConstantBasalBeta{M} <: AbstractBasalBeta
    beta::M
end

"""
$(TYPEDSIGNATURES)

A struct that defines the pseudo-plastic power law for basal sliding.
When passed to [`basal_shear_stress!`](@ref), the basal shear stress is calculated
following Eq. (24) of [robinson_description_2020](@citet). This covers plastic friction
(q = 0), linear friction (q = 1), and power-law friction (q > 1).

# Fields:
 - `v_0::T`: Reference velocity (m/yr).
 - `v_reg::T`: Regularization velocity (m/yr).
 - `q::T`: Exponent (1).
"""
@kwdef struct PseudoPlasticPowerBasalBeta{T} <: AbstractBasalBeta
    v_0::T = 1e2
    v_reg::T = 1e-3
    q::T = 0.8
end

"""
    CoulombBasalBeta{T}

A struct that defines the regularized Coulomb friction law for basal sliding.
When passed to [`basal_shear_stress!`](@ref), the basal shear stress is calculated
following Eq. (25) of [robinson_description_2020](@citet).

# Fields:
 - `v_0::T`: Reference (and regularization) velocity (m/yr).
 - `q::T`: Exponent (1).
"""
@kwdef struct CoulombBasalBeta{T} <: AbstractBasalBeta
    v_0::T = 1e2
    q::T = 0.2
end

"""
$(TYPEDSIGNATURES)
"""
function basal_beta(
    c_bed,
    v_basal,
    bb::ConstantBasalBeta,
)
    return bb.beta
end
function basal_beta(
    c_bed,
    v_basal,
    bb::PseudoPlasticPowerBasalBeta,
)
    (; v_0, v_reg, q) = bb
    v_basal_norm = norm(v_basal) + v_reg
    return c_bed * (v_basal_norm / v_0) ^ q / v_basal_norm
end
function basal_beta(
    c_bed,
    v_basal,
    bb::CoulombBasalBeta,
)
    (; v_0, q) = bb
    v_basal_norm = norm(v_basal)
    return c_bed * (v_basal_norm / (v_basal_norm + v_0)) ^ q / v_basal_norm
end

"""
$(TYPEDSIGNATURES)
"""
function basal_beta!(
    β_basal,
    c_bed,
    v_basal,
    bb::AbstractBasalBeta,
)
    map!((c, v) -> basal_beta(c, v, bb), β_basal, c_bed, v_basal)
    return nothing
end

#######################################################################
# Grounding zone
#######################################################################
abstract type AbstractBasalBetaGroundingZone end

"""
$(TYPEDSIGNATURES)
"""
struct FgroundBasalBetaGroundingZone <: AbstractBasalBetaGroundingZone
end

"""
$(TYPEDSIGNATURES)
"""
struct FractionBasalBetaGroundingZone <: AbstractBasalBetaGroundingZone
end

"""
$(TYPEDSIGNATURES)
"""
@kwdef struct HgroundBasalBetaGroundingZone{T} <: AbstractBasalBetaGroundingZone
    H_grounded_lim::T = 10
end

"""
$(TYPEDSIGNATURES)
"""
@kwdef struct ZstarBasalBetaGroundingZone{T} <: AbstractBasalBetaGroundingZone
    is_normalized::Bool = true
end

"""
$(TYPEDSIGNATURES)
"""
function basal_beta_gz(
    beta,
    mask_gz,
    f_grounded,
    bbgz::FgroundBasalBetaGroundingZone,
)
    @assert 0 <= f_grounded <= 1 "Grounding line fraction f_grounded must be in [0, 1]"
    if mask_gz
        return f_grounded * beta
    else
        return beta
    end
end

function basal_beta_gz(
    beta,
    mask_gz,
    f_gzone,
    bbgz::FractionBasalBetaGroundingZone,
)
    @assert 0 <= f_gzone <= 1 "Grounding zone fraction f_gzone must be in [0, 1]"
    if mask_gz
        return beta * f_gzone
    else
        return beta
    end
end

function basal_beta_gz(
    beta,
    mask_gz,
    H_grounded,
    bbgz::HgroundBasalBetaGroundingZone,
)

    if mask_gz
        f_scale = max(min(H_grounded, bbgz.H_grounded_lim) / bbgz.H_grounded_lim, 0)
    else
        f_scale = 1
    end
    return beta * f_scale
end

function basal_beta_gz(
    beta,
    H_eff,
    z_bed,
    z_sl,
    ρ_seawater_div_ρ_ice,
    bbgz::ZstarBasalBetaGroundingZone,
)
    
    if z_bed >= z_sl
        f_scale = H_eff
    else
        f_scale = max(H_eff - (z_sl - z_bed) * ρ_seawater_div_ρ_ice, 0)
    end

    if bbgz.is_normalized && H_eff > 0
        f_scale /= H_eff
    end

    return beta * f_scale
end

# TODO: apply saturation to all cases
function saturate_basal_beta(β, f_ground, β_min)
    if f_ground == 1
        return max(β, β_min)
    else
        return 0
    end
end

abstract type Topography end
abstract type Constants end

"""
$(TYPEDSIGNATURES)
"""
function basal_beta_gz!(β, topo::Topography, c::Constants, bbgz::FgroundBasalBetaGroundingZone)
    (; mask_gz, f_grounded, mask_grounded) = topo
    map!((b, m, f) -> basal_beta_gz(b, m, f, bbgz), β, mask_gz, f_grounded)
    map!((b, f, bm) -> saturate_basal_beta(b, f, bm), β, f_grounded, bbgz.β_min)
    return nothing
end
function basal_beta_gz!(β, topo::Topography, c::Constants, bbgz::FractionBasalBetaGroundingZone)
    (; mask_gz, f_gzone) = topo
    map!((b, m, f) -> basal_beta_gz(b, m, f, bbgz), β, mask_gz, f_gzone)
    return nothing
end
function basal_beta_gz!(β, topo::Topography, c::Constants, bbgz::HgroundBasalBetaGroundingZone)
    (; mask_gz, H_grounded) = topo
    map!((b, m, H) -> basal_beta_gz(b, m, H, bbgz), β, mask_gz, H_grounded)
    return nothing
end
function basal_beta_gz!(β, topo::Topography, c::Constants, bbgz::ZstarBasalBetaGroundingZone)
    (; z_bed, z_sl, H_eff) = topo
    (; ρ_seawater_div_ρ_ice) = c
    map!((b, H, z_b, z_s) -> basal_beta_gz(b, H, z_b, z_s, ρ_seawater_div_ρ_ice, bbgz),
        β, H_eff, z_bed, z_sl)
    return nothing
end

######################################################################
# Roughness sampling
#######################################################################
"""
$(TYPEDSIGNATURES)
"""
abstract type AbstractBedRoughnessSampling end

struct StddevBedRoughnessSampling{T} <: AbstractBedRoughnessSampling
    n_sigma::T
    n_sd::Int
    f_sd::Vector{T}
    w_sd::Vector{T}
    samples::Vector{T}
end

function StddevBedRoughnessSampling(T, n_sigma, n_sd)
    n_sigma = T(n_sigma)
    f_sd_min = -n_sigma
    f_sd_max = n_sigma

    f_sd = zeros(T, n_sd)
    w_sd = ones(T, n_sd)
    samples = Vector{T}(undef, n_sd)

    if n_sd > 1
        for q in 1:n_sd
            f_sd[q] = f_sd_min + (f_sd_max - f_sd_min) * (q - 1) / (n_sd - 1)
            w_sd[q] = exp(-0.5 * f_sd[q] ^ 2) / sqrt(2 * π)
        end
    end

    w_sd ./= sum(w_sd)

    return StddevBedRoughnessSampling{T}(n_sigma, n_sd, f_sd, w_sd, samples)
end

#######################################################################
# Sediment scaling
#######################################################################
abstract type AbstractSedimentScaling end

struct NoSedimentScaling <: AbstractSedimentScaling end

struct LinearSedimentScaling{T} <: AbstractSedimentScaling
    H_sed_min::T
    H_sed_max::T
end

function f_sediment(H_sediment, ss::NoSedimentScaling)
    return 0
end

# TODO: check if this is correct!
function f_sediment(H_sediment, ss::LinearSedimentScaling)
    if H_sediment < ss.H_sed_min
        return 0
    elseif H_sediment > ss.H_sed_max
        return 1
    else
        return (H_sediment - ss.H_sed_min) / (ss.H_sed_max - ss.H_sed_min)
    end
end

#######################################################################
# Reference bed friction coefficient c_bed_ref
#######################################################################
abstract type AbstractCbedRef end

struct BypassCbedRef <: AbstractCbedRef end

struct HomogeneousCbedRef{T} <: AbstractCbedRef
    c_bed_ref::T
end

struct LinearElevationCbedRef{T} <: AbstractCbedRef
    z0::T
    z1::T
end

@kwdef struct ExponentialElevationCbedRef{T} <: AbstractCbedRef
    z0::T
    z1::T
end

c_bed_ref(cf_ref, bt::BypassCbedRef) = cf_ref
c_bed_ref(cf_ref, bt::HomogeneousCbedRef) = bt.c_bed_ref

function c_bed_ref(
    cf_ref,
    z_bed,
    z_bed_σ,
    z_sl,
    H_sediment,
    brs::StddevBedRoughnessSampling,
    bt::ABT,
    ss::ASS,
) where {ABT<:AbstractCbedRef, ASS<:AbstractSedimentScaling}

    (; n_sd, f_sd, w_sd, samples) = brs
    for q in 1:n_sd
        λ_bed = lambda_bed(z_bed + f_sd[q] * z_bed_σ, z_sl, bt.z0, bt.z1, bt)
        samples[q] = max(cf_ref * λ_bed, cf_min)
    end

    return sum(samples .* w_sd) * (1 - f_sediment(H_sediment, ss))
end

# TODO: Check lambda bed. When comparing the linear and exponential case
# I feel like they are not consistent with each other.
function lambda_bed(z_bed, z_sl, z0, z1, bt::LinearElevationCbedRef)
    z_rel = z_sl - z_bed
    return saturate(
        (z_rel - z0) / (z1 - z0),
        0.0,
        1.0,
    )
end

function lambda_bed(z_bed, z_sl, z0, z1, bt::ExponentialElevationCbedRef)
    z_rel = z_sl - z_bed
    return saturate(
        1.0 - exp(- (z_rel - z0) / (z1 - z0)),
        0.0,
        1.0,
    )
end

#######################################################################
# Bed friction coefficient c_bed
#######################################################################
abstract type AbstractCbed end

struct AngularCbed <: AbstractCbed end
struct LinearCbed <: AbstractCbed end

function c_bed(c_bed_ref, N_eff, cb::AngularCbed)
    return tan( deg2rad(c_bed_ref) ) * N_eff
end

function c_bed(c_bed_ref, N_eff, cb::LinearCbed)
    return c_bed_ref * N_eff
end


#######################################################################
# legacy
#######################################################################
"""
$(TYPEDSIGNATURES)

Compute the basal shear stress `τ_basal` based on the basal velocity `v_basal`, the basal friction coefficient `c_basal`, and the basal friction model `friction<:AbstractBasalFriction`.
"""
function basal_shear_stress(
    v_basal::V,
    c_basal,
    friction::PseudoPlasticPowerBasalBeta,
) where {V<:AbstractVector{<:Real}}
    (; v_0, q) = friction
    v_basal_norm = norm(v_basal)
    return - c_basal * (v_basal_norm / v_0) ^ q * v_basal / v_basal_norm
end

function basal_shear_stress(
    v_basal::V,
    c_basal,
    friction::CoulombBasalBeta,
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
    friction::BasalFriction,
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
    friction::CoulombBasalBeta{T},
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
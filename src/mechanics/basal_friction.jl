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
abstract type AbstractFriction end

"""
$(TYPEDEF)

Concrete [`AbstractFriction`](@ref) bundling the sub-models needed to update the basal
friction via [`basal_friction!`](@ref).

# Fields:
 - `beta<:AbstractBasalBeta`: model for the basal friction coefficient.
 - `beta_gz<:AbstractBasalBetaGroundingZone`: model for the basal friction coefficient at the grounding zone.
 - `roughness_sampling<:AbstractBedRoughnessSampling`: model for the bed roughness sampling.
 - `c_bed_ref<:AbstractCbedRef`: model for the reference bed friction coefficient.
 - `c_bed<:AbstractCbed`: model for the bed friction coefficient.
"""
struct BasalFriction{
    BB,     # <: AbstractBasalBeta
    BBGZ,   # <: AbstractBasalBetaGroundingZone
    BRS,    # <: AbstractBedRoughnessSampling
    CBR,    # <: AbstractCbedRef
    CB,     # <: AbstractCbed
} <: AbstractFriction
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
    basal_friction!(τ_basal, v_basal, β_basal, c_bed, c_bed_ref, N_eff, z_bed, z_bed_σ, z_sl, bf)
    return nothing
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
    pointwise!(basal_shear_stress, τ_basal, (v_basal, β_basal))
    return nothing
end

#####################################################################
# Basal friction coefficient β
#####################################################################
"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the basal friction coefficient computation via [`basal_beta`](@ref).

!!! note "Reference velocities are `m yr⁻¹`, matching the solver's own velocities"
    The sliding laws below compare `|v_b|` against a reference velocity `v_0`, so the two
    have to share a time unit for the ratio to mean anything. Both are `m yr⁻¹` under
    Pagos's `(m, yr, Pa)` convention (see [`Constants`](@ref)), which also makes `β` come
    out in `Pa yr m⁻¹` — the unit `β F₂` needs to be dimensionless in
    [`beta_eff_diva!`](@ref).
"""
abstract type AbstractBasalBeta end

"""
$(TYPEDSIGNATURES)

Prescribed basal friction coefficient. The basal shear stress is linear in the
sliding velocity:

```math
\\begin{aligned}
\\boldsymbol{\\tau}_{\\mathrm{b}} = -\\beta \\, \\mathbf{v}_{\\mathrm{b}}
\\end{aligned}
```

# Fields
 - `beta::M`: Basal friction coefficient ``\\beta`` (scalar or array).
"""
struct PrescribedBasalBeta{M} <: AbstractBasalBeta
    beta::M
end

"""
$(TYPEDSIGNATURES)

Pseudo-plastic power-law basal sliding following Eq. (24) of [robinson_description_2020](@citet):

```math
\\begin{aligned}
\\boldsymbol{\\tau}_{\\mathrm{b}} = -c_{\\mathrm{b}}
\\left(\\frac{|\\mathbf{v}_{\\mathrm{b}}| + v_{\\mathrm{reg}}}{v_0}\\right)^{\\!q}
\\frac{\\mathbf{v}_{\\mathrm{b}}}{|\\mathbf{v}_{\\mathrm{b}}| + v_{\\mathrm{reg}}}
\\end{aligned}
```

Setting ``q = 0`` recovers perfectly plastic (Coulomb) sliding where the stress
magnitude equals ``c_{\\mathrm{b}}`` everywhere; ``q = 1`` gives linear (viscous)
sliding. Intermediate values represent sub-plastic flow typical of soft-bedded
ice streams. The regularization velocity ``v_{\\mathrm{reg}}`` prevents a singularity
at ``|\\mathbf{v}_{\\mathrm{b}}| = 0``.

# Fields
 - `v_0::T`: Reference velocity (``\\mathrm{m}\\,\\mathrm{yr}^{-1}``).
 - `v_reg::T`: Regularization velocity (``\\mathrm{m}\\,\\mathrm{yr}^{-1}``).
 - `q::T`: Sliding exponent (dimensionless).
"""
@kwdef struct PseudoPlasticPowerBasalBeta{T} <: AbstractBasalBeta
    v_0::T = 1e2
    v_reg::T = 1e-3
    q::T = 0.8
end

"""
$(TYPEDSIGNATURES)

Regularized Coulomb basal sliding following Eq. (25) of [robinson_description_2020](@citet),
with optional near-zero regularization following [zoet_slip_2020](@citet):

```math
\\begin{aligned}
\\boldsymbol{\\tau}_{\\mathrm{b}} = -c_{\\mathrm{b}}
\\left(\\frac{|\\mathbf{v}_{\\mathrm{b}}| + v_{\\mathrm{reg}}}{|\\mathbf{v}_{\\mathrm{b}}| + v_{\\mathrm{reg}} + v_0}\\right)^{\\!q}
\\frac{\\mathbf{v}_{\\mathrm{b}}}{|\\mathbf{v}_{\\mathrm{b}}| + v_{\\mathrm{reg}}}
\\end{aligned}
```

Unlike [`PseudoPlasticPowerBasalBeta`](@ref), the stress magnitude saturates to
``c_{\\mathrm{b}}`` at high sliding velocities, recovering the classical Coulomb
friction limit ``|\\boldsymbol{\\tau}_{\\mathrm{b}}| \\leq c_{\\mathrm{b}}``.
The velocity ``v_0`` sets the transition between the low-velocity growing regime
and the Coulomb plateau. Setting ``v_{\\mathrm{reg}} = 0`` (default) recovers
Eq. (25) of [robinson_description_2020](@citet); a small positive ``v_{\\mathrm{reg}}``
regularizes both the stress magnitude and the sliding direction at
``|\\mathbf{v}_{\\mathrm{b}}| = 0`` following [zoet_slip_2020](@citet).

# Fields
 - `v_0::T`: Transition velocity (``\\mathrm{m}\\,\\mathrm{yr}^{-1}``).
 - `v_reg::T`: Regularization velocity (``\\mathrm{m}\\,\\mathrm{yr}^{-1}``).
 - `q::T`: Sliding exponent (dimensionless).
"""
@kwdef struct CoulombBasalBeta{T} <: AbstractBasalBeta
    v_0::T = 1e2
    v_reg::T = 0.0
    q::T = 0.2
end

"""
$(TYPEDSIGNATURES)
"""
function basal_beta(
    c_bed,
    v_basal,
    bb::PrescribedBasalBeta,
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
    (; v_0, v_reg, q) = bb
    v_basal_norm = norm(v_basal) + v_reg
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
    pointwise!(basal_beta, β_basal, (c_bed, v_basal), (bb,))
    return nothing
end

#######################################################################
# Grounding zone
#######################################################################
abstract type AbstractBasalBetaGroundingZone end

"""
$(TYPEDSIGNATURES)
"""
@kwdef struct FgroundBasalBetaGroundingZone{T} <: AbstractBasalBetaGroundingZone
    # @dev TODO: `0.0` is a neutral placeholder (no additional floor beyond
    # `saturate_basal_beta`'s own branches) — set a real minimum once one is known.
    β_min::T = 0.0
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

# @dev TODO: apply saturation to all cases, not just `FgroundBasalBetaGroundingZone`.
function saturate_basal_beta(β, f_ground, β_min)
    if f_ground == 1
        return max(β, β_min)
    else
        return 0
    end
end


"""
$(TYPEDSIGNATURES)
"""
function basal_beta_gz!(β, topo, c::Constants, bbgz::FgroundBasalBetaGroundingZone)
    (; mask_gz, f_grounded, mask_grounded) = topo
    pointwise!(basal_beta_gz, β, (β, mask_gz, f_grounded), (bbgz,))
    pointwise!(saturate_basal_beta, β, (β, f_grounded), (bbgz.β_min,))
    return nothing
end
function basal_beta_gz!(β, topo, c::Constants, bbgz::FractionBasalBetaGroundingZone)
    (; mask_gz, f_gzone) = topo
    pointwise!(basal_beta_gz, β, (β, mask_gz, f_gzone), (bbgz,))
    return nothing
end
function basal_beta_gz!(β, topo, c::Constants, bbgz::HgroundBasalBetaGroundingZone)
    (; mask_gz, H_grounded) = topo
    pointwise!(basal_beta_gz, β, (β, mask_gz, H_grounded), (bbgz,))
    return nothing
end
function basal_beta_gz!(β, topo, c::Constants, bbgz::ZstarBasalBetaGroundingZone)
    (; z_bed, z_sl, H_eff) = topo
    (; density_seawater, density_ice) = c
    ρ_seawater_div_ρ_ice = density_seawater / density_ice
    pointwise!(basal_beta_gz, β, (β, H_eff, z_bed, z_sl), (ρ_seawater_div_ρ_ice, bbgz))
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

# @dev TODO: check if this is correct!
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
        # @dev TODO: floored at zero as a placeholder (non-negative friction) until an
        # actual minimum `cf_min` is wired in.
        samples[q] = max(cf_ref * λ_bed, zero(cf_ref))
    end

    return sum(samples .* w_sd) * (1 - f_sediment(H_sediment, ss))
end

# @dev TODO: check `lambda_bed` — the linear and exponential cases don't look consistent
# with each other.
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

"""
$(TYPEDSIGNATURES)
"""
function c_bed!(c_bed_out, c_bed_ref, N_eff, cb::AbstractCbed)
    pointwise!(c_bed, c_bed_out, (c_bed_ref, N_eff), (cb,))
    return nothing
end


#######################################################################
# legacy
#######################################################################
"""
$(TYPEDSIGNATURES)

Compute the basal shear stress `τ_basal` based on the basal velocity `v_basal`, the basal friction coefficient `c_basal`, and the basal friction model `friction<:AbstractBasalBeta`.
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
    (; v_0, v_reg, q) = friction
    v_basal_norm = norm(v_basal) + v_reg
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
    # @dev: `basal_shear_stress(::Any, ::Any, ::BasalFriction)` has no dedicated scalar
    # method — only the generic array-wrapper above and the
    # `PseudoPlasticPowerBasalBeta`/`CoulombBasalBeta` methods, neither of which matches a
    # scalar call with a `BasalFriction`. Untested (nothing in test/ exercises
    # basal_friction.jl); possibly already broken.
    pointwise!(basal_shear_stress, τ_basal, (v_basal, c_basal), (friction,))
    return nothing
end

"""
$(TYPEDSIGNATURES)

This will implement Eq. 2 of Zoet & Iverson (2020).
"""
function ut_with_clast()
    error("Not implemented yet")
end
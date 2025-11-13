###########################################################
# Structs
###########################################################

"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the effective pressure computation via [`effective_pressure`](@ref).
"""
abstract type AbstractEffectivePressure{T<:AbstractFloat} end

"""
$(TYPEDSIGNATURES)

A struct that passes an externally computed field effective pressure field via [`effective_pressure`](@ref).

# Fields
 - `p::Matrix{T}`: effective pressure field (Pa).
"""
struct ExternalEffectivePressure{T} <: AbstractEffectivePressure{T}
    p::Matrix{T}
end

"""
$(TYPEDSIGNATURES)

A struct that passes a constant effective pressure value via [`effective_pressure`](@ref).

# Fields
 - `p::T`: effective pressure value (Pa).
"""
struct ConstantEffectivePressure{T} <: AbstractEffectivePressure{T}
    p::T
end

"""
$(TYPEDSIGNATURES)

A struct that computes the overburden pressure as effective pressure via [`effective_pressure`](@ref).

# Fields
 - `ρ::T`: ice density (``\\mathrm{kg \\, m^{-3}}``).
"""
@kwdef struct OverburdenEffectivePressure{T} <: AbstractEffectivePressure{T}
    ρ::T = 918.0      # ice density
    g::T = 9.81       # gravitational acceleration
end

"""
$(TYPEDSIGNATURES)

A struct that computes the Leguy et al. (2020) effective pressure via [`effective_pressure`](@ref).

# Fields
 - `p::T`: marine connectivity exponent (0: none, 1: full).
"""
struct LeguyEffectivePressure{T} <: AbstractEffectivePressure{T}
    p::T
    overburden::OverburdenEffectivePressure{T}
end

"""
$(TYPEDSIGNATURES)

A struct that computes the till pressure via [`effective_pressure`](@ref), following Bueler and van Pelt (2015).

# Fields
 - `H_w_max::T=2.0`: saturation water thickness (m).
 - `N0::T=1e3`: reference effective pressure (Pa).
 - `delta::T=0.04`: fraction of overburden pressure for saturated till.
 - `e0::T=0.69`: reference void ratio at N0.
 - `Cc::T=0.12`: till compressibility.
"""
@kwdef struct TillEffectivePressure{T} <: AbstractEffectivePressure{T}
    H_w_max::T = 2.0        # saturation water thickness
    N0::T = 1e3             # reference effective pressure
    delta::T = 0.04         # fraction of overburden pressure for saturated till
    e0::T = 0.69            # reference void ratio at N0
    Cc::T = 0.12            # till compressibility
    overburden::OverburdenEffectivePressure{T} = OverburdenEffectivePressure{T}()
end

###########################################################
# Dispatch
###########################################################

"""
$(TYPEDSIGNATURES)

Compute the effective pressure using the specified effective pressure model `eff_pressure<:AbstractEffectivePressure`.
"""
function effective_pressure(
    H_ice,
    H_water,
    eff_pressure::ExternalEffectivePressure,
)
    return eff_pressure.p
end

function effective_pressure(
    H_ice::R1,
    H_water::R2,
    eff_pressure::ConstantEffectivePressure,
) where {R1<:Real, R2<:Real}
    return eff_pressure.p
end

function effective_pressure(
    H_ice::R1,
    H_water::R2,
    eff_pressure::OverburdenEffectivePressure,
) where {R1<:Real, R2<:Real}
    return eff_pressure.ρ * eff_pressure.g * H_ice
end

function effective_pressure(
    H_ice,
    H_water,
    eff_pressure::LeguyEffectivePressure,
)
    error("Not implemented yet")
end

function effective_pressure(
    H_ice,
    H_water,
    eff_pressure::TillEffectivePressure,
)
    error("Not implemented yet")
end
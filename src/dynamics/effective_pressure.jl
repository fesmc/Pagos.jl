###########################################################
# Structs
###########################################################
"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the effective pressure computation via [`effective_pressure`](@ref).
"""
abstract type AbstractEffectivePressure end

"""
$(TYPEDSIGNATURES)

A struct that passes a constant effective pressure value via [`effective_pressure`](@ref).

# Fields
 - `p::M`: effective pressure value (Pa). Can be a scalar or an array.
"""
struct PrescribedEffectivePressure{M} <: AbstractEffectivePressure
    p::M
end

"""
$(TYPEDSIGNATURES)

A struct that computes the overburden pressure as effective pressure via [`effective_pressure`](@ref).

# Fields
 - `ρ::T`: ice density (``\\mathrm{kg \\, m^{-3}}``).
"""
@kwdef struct OverburdenEffectivePressure{T} <: AbstractEffectivePressure
    ρ::T = 918.0      # ice density
    g::T = 9.81       # gravitational acceleration
end

"""
$(TYPEDSIGNATURES)

A struct that computes the Leguy et al. (2020) effective pressure via [`effective_pressure`](@ref).

# Fields
 - `p::T`: marine connectivity exponent (0: none, 1: full).
"""
struct LeguyEffectivePressure{T} <: AbstractEffectivePressure
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
@kwdef struct TillEffectivePressure{T} <: AbstractEffectivePressure
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
    H_eff,
    f_ground,
    eff_pressure::AbstractEffectivePressure,
)
    if f_ground > 0.0
        return effective_pressure(H_eff, eff_pressure)
    else
        return 0
    end
end
function effective_pressure(
    H_eff,
    eff_pressure::PrescribedEffectivePressure,
)
    return eff_pressure.p
end
function effective_pressure(
    H_eff,
    eff_pressure::OverburdenEffectivePressure,
)
    return eff_pressure.ρ * eff_pressure.g * H_eff
end
function effective_pressure(
    H_eff,
    eff_pressure::LeguyEffectivePressure,
)
    error("Not implemented yet")
end
function effective_pressure(
    H_eff,
    eff_pressure::TillEffectivePressure,
)
    error("Not implemented yet")
end

"""
$(TYPEDSIGNATURES)

Same as [`effective_pressure`](@ref) but operates in place.
"""
function effective_pressure!(N_eff, H_eff, eff_pressure::PrescribedEffectivePressure)
    return nothing
end
function effective_pressure!(N_eff, H_eff, eff_pressure)
    map!((h_eff) -> effective_pressure(h_eff, eff_pressure), N_eff, H_eff)
end

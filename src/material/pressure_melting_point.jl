"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the pressure melting point computation via [`pressure_melting_point`](@ref), [`pressure_melting_point!`](@ref), [`relative_temperature`](@ref) and [`relative_temperature!`](@ref).
"""
abstract type AbstractPressureMeltingPoint end

"""
$(TYPEDSIGNATURES)

Prescribed pressure melting point, independent of pressure. Useful for testing, and for cold-ice models where the pressure correction is negligible.

# Fields
 - `T_0::T=273.15`: melting point (``\\mathrm{K}``).
"""
@kwdef struct PrescribedPressureMeltingPoint{T} <: AbstractPressureMeltingPoint
    T_0::T = 273.15
end

"""
$(TYPEDSIGNATURES)

A linear pressure melting point parameterization:

```math
\\begin{align}
T_{\\mathrm{m}} = T_{0} - \\beta \\, p
\\end{align}
```

where ``T_{0}`` is the melting point at standard pressure (273.15 K), ``\\beta`` is the Clausius-Clapeyron constant, and ``p`` is the pressure. The temperature relative to the pressure melting point is then given by:

```math
\\begin{align}
T' = T + \\beta \\, p
\\end{align}
```

For pure ice, the Clausius-Clapeyron constant yields ``\\beta = 7.42 \\times 10^{-8} \\, \\mathrm{K}\\,\\mathrm{Pa}^{-1}``, but under realistic conditions, the value for air-saturated ice is closer to ``\\beta = 9.8 \\times 10^{-8} \\, \\mathrm{K}\\,\\mathrm{Pa}^{-1}`` ([greve_dynamics_2009](@citet), p. 53-54). Therefore, we use the latter as default value.

# Fields
 - `T_0::T=273.15`: melting point at standard pressure (``\\mathrm{K}``).
 - `β::T=9.8e-8`: Clausius-Clapeyron constant (``\\mathrm{K}\\,\\mathrm{Pa}^{-1}``).
"""
@kwdef struct LinearPressureMeltingPoint{T} <: AbstractPressureMeltingPoint
    T_0::T = 273.15     # K, melting point at standard pressure
    β::T = 9.8e-8       # K Pa^-1, Greve and Blatter (2009), p. 54, (Hooke 2005)
end

"""
$(TYPEDSIGNATURES)

Linearized pressure- and salinity-dependent melting point at the ice–ocean interface, following [holland_thermohaline_1999](@citet) and the ISOMIP+ protocol
[asay-davis_experimental_2016](@citet):

```math
\\begin{align}
T_{\\mathrm{f}} = \\lambda_1 \\, S + \\lambda_2 + \\lambda_3 \\, p
\\end{align}
```

where ``S`` is salinity (psu), ``p`` is pressure (Pa), and the default coefficients are those of the ISOMIP+ protocol. Use this parameterization at the ice–ocean interface (ice shelf base, grounding line) wherever a melt-rate parameterization (PICO, PICOP, plume model) requires the local freezing point.

# Fields
 - `λ₁::T=-0.0573`: salinity coefficient (``\\mathrm{K}\\,\\mathrm{psu}^{-1}``).
 - `λ₂::T=0.0832`: constant offset (``\\mathrm{K}``).
 - `λ₃::T=-7.53e-8`: pressure coefficient (``\\mathrm{K}\\,\\mathrm{Pa}^{-1}``).
"""
@kwdef struct LinearSalinityPressureMeltingPoint{T} <: AbstractPressureMeltingPoint
    λ₁::T = -0.0573     # K psu⁻¹, salinity coefficient (ISOMIP+)
    λ₂::T = 0.0832      # K, constant offset (Holland & Jenkins 1999)
    λ₃::T = -7.53e-8    # K Pa⁻¹, pressure coefficient (ISOMIP+)
end

"""
$(TYPEDSIGNATURES)

Get the melting point `T_m` at pressure `p` based on the pressure melting point parameterization `law<:AbstractPressureMeltingPoint`.
"""
function pressure_melting_point(p, law::PrescribedPressureMeltingPoint)
    return law.T_0
end

function pressure_melting_point(p, law::LinearPressureMeltingPoint)
    return law.T_0 - law.β * p
end

function pressure_melting_point(p, S, law::LinearSalinityPressureMeltingPoint)
    return law.λ₂ + law.λ₁ * S + law.λ₃ * p
end

"""
$(TYPEDSIGNATURES)

Same as [`pressure_melting_point`](@ref) but operates in place.
"""
function pressure_melting_point!(Tm, p, law)
    map!(x -> pressure_melting_point(x, law), Tm, p)
    return
end

function pressure_melting_point!(Tf, p, S, law::LinearSalinityPressureMeltingPoint)
    map!((pp, ss) -> pressure_melting_point(pp, ss, law), Tf, p, S)
    return
end

"""
$(TYPEDSIGNATURES)

Get the relative temperature `T'` based on the absolute temperature `T` and pressure `p` using the pressure melting point parameterization `law<:AbstractPressureMeltingPoint`.
"""
function relative_temperature(T, p, law::PrescribedPressureMeltingPoint)
    return T
end

function relative_temperature(T, p, law::LinearPressureMeltingPoint)
    return T + law.β * p
end

"""
$(TYPEDSIGNATURES)

Same as [`relative_temperature`](@ref) but operates in place.
"""
function relative_temperature!(Tprime, T, p, law::LinearPressureMeltingPoint)
    map!(x -> relative_temperature(x[1], x[2], law), Tprime, T, p)
    return
end

"""
$(TYPEDSIGNATURES)

Get the thermal forcing ``\\Delta T = T - T_f(S, p)`` at the ice–ocean interface.
This is the ocean analogue of [`relative_temperature`](@ref) and is the primary driver
of basal melting in PICO, PICOP, and plume-type parameterizations.
"""
function thermal_forcing(T, p, S, law::LinearSalinityPressureMeltingPoint)
    return T - pressure_melting_point(p, S, law)
end

"""
$(TYPEDSIGNATURES)

Same as [`thermal_forcing`](@ref) but operates in place.
"""
function thermal_forcing!(ΔT, T, p, S, law::LinearSalinityPressureMeltingPoint)
    map!((ti, pi, si) -> thermal_forcing(ti, pi, si, law), ΔT, T, p, S)
    return
end

# function update_relative_temperature!(
#     T_rel_ice::Array{T, 3},
#     T_ice::Array{T, 3},
#     z::Vector{T},
#     H::Matrix{T},
#     pressure_melting_point::LinearPressureMeltingPoint{T},
#     c::PhysicalConstants{T},
#     idx::Matrix{CartesianIndex{2}},
# ) where {T<:AbstractFloat}
    
#     β = pressure_melting_point.β
#     (; ρ_ice, g) = c
#     for i in idx
#         for l in axes(T_rel_ice, 3)
#             T_rel_ice[i, l] = T_ice[i, l] + β * pressure_column(ρ_ice, g, z[l] * H[i, l])
#         end
#     end

#     return nothing
# end
"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the pressure melting point computation
via [`melting_point`](@ref), [`update_melting_point`](@ref), [`get_relative_temperature`](@ref) and
[`update_relative_temperature!`](@ref).
"""
abstract type AbstractPressureMeltingPoint end

"""
$(TYPEDSIGNATURES)

A linear pressure melting point parameterization:

```math
\\begin{align}
T_{m} = T_{0} - \\beta p
\\end{align}
```

where ``T_{0}`` is the melting point at standard pressure (273.15 K), ``\\beta`` is the Clausius-Clapeyron constant, and ``p`` is the pressure. The temperature relative to the pressure melting point is then given by:

```math
\\begin{align}
T' = T + \\beta p
\\end{align}
```

For pure ice, the Clausius-Clapeyron constant yields ``\\beta = 7.42 \\times 10^{-8} \\, \\mathrm{K Pa^{-1}}``, but under realistic conditions, the value for air-saturated ice is closer to ``\\beta = 9.8 \\times 10^{-8} \\, \\mathrm{K Pa^{-1}}`` ([greve_dynamics_2009](@citet), p. 53-54). Therefore, we use the latter as default value.

# Fields
 - `T_0::T=273.15`: melting point at standard pressure (``\\mathrm{K}``).
 - `β::T=9.8e-8`: Clausius-Clapeyron constant (``\\mathrm{K \\, Pa^{-1}}``).
"""
@kwdef struct LinearPressureMeltingPoint{T} <: AbstractPressureMeltingPoint
    T_0::T = 273.15     # K, melting point at standard pressure
    β::T = 9.8e-8       # K Pa^-1, Greve and Blatter (2009), p. 54, (Hooke 2005)
end

"""
$(TYPEDSIGNATURES)

Get the melting point `T_m` at pressure `p` based on the pressure melting point parameterization `pmp<:AbstractPressureMeltingPoint`.
"""
function pressure_melting_point(p, pmp::LinearPressureMeltingPoint)
    return pmp.T_0 - pmp.β * p
end

"""
$(TYPEDSIGNATURES)

Same as [`pressure_melting_point`](@ref) but operates in place.
"""
function pressure_melting_point!(Tm, p, pmp::LinearPressureMeltingPoint)
    map!(x -> pressure_melting_point(x, pmp), Tm, p)
    return
end

"""
$(TYPEDSIGNATURES)

Get the relative temperature `T'` based on the absolute temperature `T` and pressure `p` using the pressure melting point parameterization `pmp<:AbstractPressureMeltingPoint`.
"""
function relative_temperature(T, p, pmp::LinearPressureMeltingPoint)
    return T + pmp.β * p
end

"""
$(TYPEDSIGNATURES)

Update the relative temperature `T'` based on the absolute temperature `T` and pressure `p` using the pressure melting point parameterization `pmp<:AbstractPressureMeltingPoint`.
"""
function relative_temperature!(Tprime, T, p, pmp::LinearPressureMeltingPoint)
    map!(x -> relative_temperature(x[1], x[2], pmp), Tprime, T, p)
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
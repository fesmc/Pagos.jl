###########################################################
# Structs
###########################################################

"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the flow law computation via [`viscosity`](@ref) and [`viscosity!`](@ref).
"""
abstract type AbstractFlowLaw end

"""
$(TYPEDSIGNATURES)

Flow law with constant viscosity.
# Fields
 - `η::T`: constant viscosity.
"""
struct ConstantViscosityFlowLaw{T} <: AbstractFlowLaw
    η::T
end

"""
$(TYPEDSIGNATURES)

Flow law that combines an [`AbstractRateFactor`](@ref) and an [`AbstractCreep`](@ref)
to compute the ice viscosity:

```math
\\begin{aligned}
\\eta(T', \\sigma_e) = \\frac{1}{2 \\, A(T') \\, f(\\sigma_e)}
\\end{aligned}
```

where ``A(T')`` is the rate factor and ``f(\\sigma_e)`` is the creep function.

# Fields
 - `rate_factor::RF`: rate factor model (``<:`` [`AbstractRateFactor`](@ref)).
 - `creep::C`: creep function model (``<:`` [`AbstractCreep`](@ref)).
"""
struct RateCreepFlowLaw{
    RF,     # <: AbstractRateFactor,
    C,     # <: AbstractCreep,
} <: AbstractFlowLaw
    rate_factor::RF
    creep::C
end

"""
$(TYPEDSIGNATURES)

Convenience function to create a Glen-Nye flow law with [`ArrheniusRateFactor`](@ref) and [`GlenNyeCreep`](@ref).
"""
function GlenNyeFlowLaw()
    rf = ArrheniusRateFactor()
    c = GlenNyeCreep()
    return RateCreepFlowLaw(rf, c)
end

"""
$(TYPEDSIGNATURES)

Convenience function to create a Regularized Glen-Nye flow law with [`ArrheniusRateFactor`](@ref) and [`RegularizedGlenNyeCreep`](@ref).
"""
function RegularizedGlenNyeFlowLaw()
    rf = ArrheniusRateFactor()
    c = RegularizedGlenNyeCreep()
    return RateCreepFlowLaw(rf, c)
end

"""
$(TYPEDSIGNATURES)

Convenience function to create a Smith-Morland flow law with [`SmithMorlandRateFactor`](@ref) and [`SmithMorlandCreep`](@ref).
"""
function SmithMorlandFlowLaw()
    rf = SmithMorlandRateFactor()
    c = SmithMorlandCreep()
    return RateCreepFlowLaw(rf, c)
end

###########################################################
# Dispatch
###########################################################

"""
$(TYPEDSIGNATURES)

Get the viscosity `η` based on the rate factor `A`, creep function `f`, and flow law parameterization `law<:AbstractFlowLaw`.
"""
function viscosity(
    A,
    f,
    law::ConstantViscosityFlowLaw,
)
    return law.η
end

function viscosity(
    A,  # rate factor
    f,  # creep function
    law::RateCreepFlowLaw,
)
    return 0.5 ./ (A .* f)
end

function viscosity(
    A::M,
    f::M,
    law,
) where {M<:AbstractArray}
    η = similar(A)
    viscosity!(η, A, f, law)
    return η
end

"""
$(TYPEDSIGNATURES)

Same as [`viscosity`](@ref) but operates in place.
"""
function viscosity!(η, A, f, law)
    map!((a, s) -> viscosity(a, s, law), η, A, f)
    return nothing
end
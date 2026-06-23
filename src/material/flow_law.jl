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
struct PrescribedViscosityFlowLaw{T} <: AbstractFlowLaw
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

Convenience function to create a Smith-Morland flow law with [`SmithMorlandRateFactor`](@ref) and [`SmithMorlandCreep`](@ref).
"""
function SmithMorlandFlowLaw()
    rf = SmithMorlandRateFactor()
    c = SmithMorlandCreep()
    return RateCreepFlowLaw(rf, c)
end

"""
$(TYPEDSIGNATURES)

Convenience function to create the [fan_flow_2025](@citet) low-strain (1–2%)
three-component flow law for a given grain size `d` (m) and absolute temperature
`T` (K). The Arrhenius pre-factors for each component are evaluated at `T` via
[`FanLowStrainGSIRateFactor`](@ref), [`FanLowStrainGSS1RateFactor`](@ref), and
[`FanLowStrainGSS2RateFactor`](@ref) and stored in a [`FanLowStrainCreep`](@ref)
struct. Because all temperature dependence is pre-baked into the creep struct, the
returned [`RateCreepFlowLaw`](@ref) uses a unit rate factor ``A = 1``.

This law is applicable to isotropic ice at low strain, or to anisotropic ice
when the deformation kinematics differ from those that formed the crystallographic
preferred orientation (for example, borehole closure, grounding-line flexure).
"""
function FanLowStrainFlowLaw(d, T)
    A_GSI  = rate_factor(T, FanLowStrainGSIRateFactor())
    A_GSS1 = rate_factor(T, FanLowStrainGSS1RateFactor())
    A_GSS2 = rate_factor(T, FanLowStrainGSS2RateFactor())
    c = FanLowStrainCreep(; d, A_GSI, A_GSS1, A_GSS2)
    return RateCreepFlowLaw(PrescribedRateFactor(one(typeof(A_GSI))), c)
end

"""
$(TYPEDSIGNATURES)

Convenience function to create the [fan_flow_2025](@citet) high-strain (≥8%)
one-component GSI flow law with [`FanHighStrainGSIRateFactor`](@ref) and
[`GlenNyeCreep`](@ref) (``n = 3.5``). Suitable for large-scale ice-sheet
modelling where ice has reached microstructural steady state (tertiary creep /
flow stress). Use ``n = 3.5`` when modelling includes temperatures above ``-5``°C;
for colder domains ``n \\approx 4`` is a better fit.
"""
function FanHighStrainFlowLaw()
    rf = FanHighStrainGSIRateFactor()
    c = GlenNyeCreep(n = 3.5)
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
    _,
    _,
    law::PrescribedViscosityFlowLaw,
)
    return law.η
end

function viscosity(
    A,
    f,
    ::RateCreepFlowLaw,
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
function viscosity!(η::AbstractVector, A, f, law)
    @tullio η[i] = viscosity(A[i], f[i], law)
    return nothing
end
function viscosity!(η::AbstractMatrix, A, f, law)
    @tullio η[i, j] = viscosity(A[i, j], f[i, j], law)
    return nothing
end
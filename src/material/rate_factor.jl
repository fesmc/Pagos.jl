###########################################################
# Structs
###########################################################

"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the rate factor computation via [`rate_factor`](@ref).
"""
abstract type AbstractRateFactor end

"""
$(TYPEDSIGNATURES)

Rate factor for ice viscosity following a constant value.

# Fields
- `A::T=1e-16`: rate factor.
"""
@kwdef struct ConstantRateFactor{T} <: AbstractRateFactor
    A::T = 1e-16
end

"""
$(TYPEDSIGNATURES)

Rate factor for ice viscosity following the piecewise Arrhenius law of [paterson_physics_1994](@citet):

```math
\\begin{aligned}
A(T') =
\\begin{cases}
E_f \\, A_{0,1} \\, e^{-Q_1 / (R \\, T')} & T' \\leq T^* \\\\
E_f \\, A_{0,2} \\, e^{-Q_2 / (R \\, T')} & T' > T^*
\\end{cases}
\\end{aligned}
```

where ``T'`` is the temperature relative to the pressure melting point (see
[`LinearPressureMeltingPoint`](@ref)) and ``T^*`` is the breakpoint temperature
below which the lower activation energy ``Q_1`` applies.

# Fields
 - `E_f::T=1.0`: enhancement factor.
 - `T_p1_p2::T=263.15`: breakpoint temperature ``T^*`` (``\\mathrm{K}``).
 - `A_0_p1::T=3.985e-13`: pre-exponential factor ``A_{0,1}`` (``\\mathrm{Pa}^{-3}\\,\\mathrm{s}^{-1}``).
 - `A_0_p2::T=1.916e3`: pre-exponential factor ``A_{0,2}`` (``\\mathrm{Pa}^{-3}\\,\\mathrm{s}^{-1}``).
 - `Q_a_p1::T=60e3`: activation energy ``Q_1`` (``\\mathrm{J}\\,\\mathrm{mol}^{-1}``).
 - `Q_a_p2::T=139e3`: activation energy ``Q_2`` (``\\mathrm{J}\\,\\mathrm{mol}^{-1}``).
 - `R::T=8.314`: universal gas constant (``\\mathrm{J}\\,\\mathrm{K}^{-1}\\,\\mathrm{mol}^{-1}``).

# Examples
```jldoctest
arrhenius_rate_factor = ArrheniusRateFactor()
T_relative_kelvin = range(-50, stop = 10, step = 0.1) .+ 273.15
A = rate_factor(T_relative_kelvin, arrhenius_rate_factor)
fig = plot_rate_factor(T_relative_kelvin .- 273.15, A)
```
"""
@kwdef struct ArrheniusRateFactor{T} <: AbstractRateFactor
    E_f::T = 1.0
    T_p1_p2::T = 263.15     # breakpoint temperature following Patterson (1994)
    A_0_p1::T = 3.985e-13   # "s^-1 Pa^-3" pre-exponential factor
    A_0_p2::T = 1.916e3     # piecewise definition following Patterson (1994)
    Q_a_p1::T = 60e3        # "J mol^-1" activation energy
    Q_a_p2::T = 139e3       #  piecewise definition following Patterson (1994)
    R::T = 8.314            # "J K^-1 mol^-1" universal gas constant
end

"""
$(TYPEDSIGNATURES)

Warning: This is not fully supported yet, since the input fields are different than those used in Glen's flow law, which is much more common.

Rate factor for ice viscosity following [smith_viscous_1981](@citet).

# Fields
 - `p1::T=0.7242`
 - `p2::T=0.3438`
 - `e1::T=11.9567`
 - `e2::T=2.9494`
"""
@kwdef struct SmithMorlandRateFactor{T} <: AbstractRateFactor
    p1::T = 0.7242
    p2::T = 0.3438
    e1::T = 11.9567
    e2::T = 2.9494
    T_0::T = 273.15
    ΔT::T = 20.0
end

"""
$(TYPEDSIGNATURES)

Rate factor for ice viscosity following [hooke_flow_1981](@citet):

```math
\\begin{aligned}
A(T') = E_f \\, A_0 \\, \\exp\\!\\left( -\\frac{Q_a}{R \\, T'} + \\frac{3C}{(T_r - T')^k} \\right)
\\end{aligned}
```

A single, continuous Arrhenius expression augmented by the proximity-to-melting
term ``3C/(T_r - T')^k``, which produces a steep upturn as ``T'`` approaches the
pressure melting point ``T_r``. Avoids the discontinuity of the piecewise
[`ArrheniusRateFactor`](@ref). Valid for ``T' < T_r``.

# Fields
 - `E_f::T=1.0`: enhancement factor.
 - `A_0::T=9.302e-7`: pre-exponential factor (``\\mathrm{Pa}^{-3}\\,\\mathrm{s}^{-1}``).
 - `Q_a::T=78.8e3`: activation energy (``\\mathrm{J}\\,\\mathrm{mol}^{-1}``).
 - `C::T=0.16612`: proximity-to-melting coefficient (``\\mathrm{K}^{k}``).
 - `T_r::T=273.39`: reference melting temperature (``\\mathrm{K}``).
 - `k::T=1.17`: proximity-to-melting exponent.
 - `R::T=8.314`: universal gas constant (``\\mathrm{J}\\,\\mathrm{K}^{-1}\\,\\mathrm{mol}^{-1}``).
"""
@kwdef struct HookeRateFactor{T} <: AbstractRateFactor
    E_f::T = 1.0
    A_0::T = 9.302e-7
    Q_a::T = 78.8e3
    C::T = 0.16612
    T_r::T = 273.39
    k::T = 1.17
    R::T = 8.314
end

"""
$(TYPEDSIGNATURES)

Rate factor for temperate ice following [lliboutry_various_1985](@citet):

```math
\\begin{aligned}
A(T', \\omega) = A_{\\mathrm{cold}}(T') \\, (1 + \\gamma \\, \\omega)
\\end{aligned}
```

where ``A_{\\mathrm{cold}}(T')`` is the piecewise [`ArrheniusRateFactor`](@ref)
evaluated at the same breakpoint temperature, and ``\\omega`` is the volumetric
liquid-water content. At ``\\omega = 0`` the law is identical to
[`ArrheniusRateFactor`](@ref); increasing ``\\omega`` shifts the entire rate-factor
curve upward proportionally. Assumes spatially uniform water content; for
spatially varying fields, call the scalar method pointwise.

# Fields
 - `ω::T=0.0`: volumetric water content fraction (dimensionless, ``0 \\leq \\omega \\leq 1``).
 - `γ::T=181.25`: water-content enhancement coefficient following [lliboutry_various_1985](@citet).
 - `E_f::T=1.0`: additional enhancement factor.
 - `T_p1_p2::T=263.15`: breakpoint temperature (``\\mathrm{K}``).
 - `A_0_p1::T=3.985e-13`: pre-exponential factor below breakpoint (``\\mathrm{Pa}^{-3}\\,\\mathrm{s}^{-1}``).
 - `A_0_p2::T=1.916e3`: pre-exponential factor above breakpoint.
 - `Q_a_p1::T=60e3`: activation energy below breakpoint (``\\mathrm{J}\\,\\mathrm{mol}^{-1}``).
 - `Q_a_p2::T=139e3`: activation energy above breakpoint.
 - `R::T=8.314`: universal gas constant (``\\mathrm{J}\\,\\mathrm{K}^{-1}\\,\\mathrm{mol}^{-1}``).
"""
@kwdef struct LliboutryDuvalRateFactor{T} <: AbstractRateFactor
    ω::T = 0.0
    γ::T = 181.25
    E_f::T = 1.0
    T_p1_p2::T = 263.15
    A_0_p1::T = 3.985e-13
    A_0_p2::T = 1.916e3
    Q_a_p1::T = 60e3
    Q_a_p2::T = 139e3
    R::T = 8.314
end

"""
$(TYPEDSIGNATURES)

Arrhenius rate factor for the grain-size insensitive (GSI) dislocation-creep
component of the low-strain (1–2%) multicomponent flow law of
[fan_flow_2025](@citet):

```math
\\begin{aligned}
A_{\\mathrm{GSI}}(T) = A_0 \\, \\exp\\!\\left(-\\frac{Q}{R \\, T}\\right)
\\end{aligned}
```

Default parameter values are the posterior medians from Bayesian inference over
305 low-strain data points (Table 1 of [fan_flow_2025](@citet)), converted from
MPa to Pa via ``\\log_{10} A_{\\mathrm{Pa}} = \\log_{10} A_{\\mathrm{MPa}} - 6n``.

# Fields
 - `A_0::T=2.951e-17`: pre-exponential factor (``\\mathrm{Pa}^{-3.6}\\,\\mathrm{s}^{-1}``).
 - `Q::T=62e3`: activation energy (``\\mathrm{J}\\,\\mathrm{mol}^{-1}``).
 - `R::T=8.314`: universal gas constant (``\\mathrm{J}\\,\\mathrm{K}^{-1}\\,\\mathrm{mol}^{-1}``).
"""
@kwdef struct FanLowStrainGSIRateFactor{T} <: AbstractRateFactor
    A_0::T = 2.951e-17   # 10^(5.07 - 6*3.6) Pa^-3.6 s^-1
    Q::T = 62e3
    R::T = 8.314
end

"""
$(TYPEDSIGNATURES)

Arrhenius rate factor for the first grain-size sensitive (GSS1) disGBS component
of the low-strain (1–2%) multicomponent flow law of [fan_flow_2025](@citet):

```math
\\begin{aligned}
A_{\\mathrm{GSS1}}(T) = A_0 \\, \\exp\\!\\left(-\\frac{Q}{R \\, T}\\right)
\\end{aligned}
```

Default parameter values are the posterior medians from Bayesian inference over
305 low-strain data points (Table 1 of [fan_flow_2025](@citet)), converted from
MPa to Pa.

# Fields
 - `A_0::T=4.677e-13`: pre-exponential factor (``\\mathrm{Pa}^{-1.9}\\,\\mathrm{m}^{1.2}\\,\\mathrm{s}^{-1}``).
 - `Q::T=52e3`: activation energy (``\\mathrm{J}\\,\\mathrm{mol}^{-1}``).
 - `R::T=8.314`: universal gas constant (``\\mathrm{J}\\,\\mathrm{K}^{-1}\\,\\mathrm{mol}^{-1}``).
"""
@kwdef struct FanLowStrainGSS1RateFactor{T} <: AbstractRateFactor
    A_0::T = 4.677e-13   # 10^(-0.93 - 6*1.9) Pa^-1.9 m^1.2 s^-1
    Q::T = 52e3
    R::T = 8.314
end

"""
$(TYPEDSIGNATURES)

Arrhenius rate factor for the second grain-size sensitive (GSS2) disGBS component
of the low-strain (1–2%) multicomponent flow law of [fan_flow_2025](@citet):

```math
\\begin{aligned}
A_{\\mathrm{GSS2}}(T) = A_0 \\, \\exp\\!\\left(-\\frac{Q}{R \\, T}\\right)
\\end{aligned}
```

The high activation energy ``Q = 182`` kJ mol⁻¹ causes this component to dominate
near the pressure melting point, capturing the steep increase in apparent ``Q``
observed experimentally above approximately ``-10`` °C without imposing a
discontinuous temperature threshold. Default parameter values are the posterior
medians from Bayesian inference over 305 low-strain data points (Table 1 of
[fan_flow_2025](@citet)), converted from MPa to Pa.

# Fields
 - `A_0::T=4.571e7`: pre-exponential factor (``\\mathrm{Pa}^{-2.5}\\,\\mathrm{m}^{1.9}\\,\\mathrm{s}^{-1}``).
 - `Q::T=182e3`: activation energy (``\\mathrm{J}\\,\\mathrm{mol}^{-1}``).
 - `R::T=8.314`: universal gas constant (``\\mathrm{J}\\,\\mathrm{K}^{-1}\\,\\mathrm{mol}^{-1}``).
"""
@kwdef struct FanLowStrainGSS2RateFactor{T} <: AbstractRateFactor
    A_0::T = 4.571e7     # 10^(22.66 - 6*2.5) Pa^-2.5 m^1.9 s^-1
    Q::T = 182e3
    R::T = 8.314
end

"""
$(TYPEDSIGNATURES)

Arrhenius rate factor for the high-strain (≥8%) one-component GSI flow law of
[fan_flow_2025](@citet):

```math
\\begin{aligned}
A_{\\mathrm{GSI}}(T) = A_0 \\, \\exp\\!\\left(-\\frac{Q}{R \\, T}\\right)
\\end{aligned}
```

Calibrated to tertiary-creep / flow-stress data (160 data points). Stress
exponent ``n = 3.5`` (or ``n \\approx 4`` when ``T < -5`` °C) gives a better
fit to the high-strain experimental record than the classic Glen exponent.
Default parameter values are the posterior medians from Table 1 of
[fan_flow_2025](@citet), converted from MPa to Pa.

# Fields
 - `A_0::T=7.943e-10`: pre-exponential factor (``\\mathrm{Pa}^{-3.5}\\,\\mathrm{s}^{-1}``).
 - `Q::T=90e3`: activation energy (``\\mathrm{J}\\,\\mathrm{mol}^{-1}``).
 - `R::T=8.314`: universal gas constant (``\\mathrm{J}\\,\\mathrm{K}^{-1}\\,\\mathrm{mol}^{-1}``).
"""
@kwdef struct FanHighStrainGSIRateFactor{T} <: AbstractRateFactor
    A_0::T = 7.943e-10   # 10^(11.9 - 6*3.5) Pa^-3.5 s^-1
    Q::T = 90e3
    R::T = 8.314
end

###########################################################
# Functions
###########################################################

"""
$(TYPEDSIGNATURES)

Get the rate factor `A` based on the temperature relative to the pressure melt point `T_relative` and the rate factor parameterization `arf<:AbstractRateFactor`.
"""
function rate_factor(
    ::T,
    crf::ConstantRateFactor,
) where {T<:Real}
    return crf.A
end

function rate_factor(
    T_relative::T,
    rf::FanLowStrainGSIRateFactor,
) where {T<:Real}
    (; A_0, Q, R) = rf
    return A_0 * exp(-Q / (R * T_relative))
end

function rate_factor(
    T_relative::T,
    rf::FanLowStrainGSS1RateFactor,
) where {T<:Real}
    (; A_0, Q, R) = rf
    return A_0 * exp(-Q / (R * T_relative))
end

function rate_factor(
    T_relative::T,
    rf::FanLowStrainGSS2RateFactor,
) where {T<:Real}
    (; A_0, Q, R) = rf
    return A_0 * exp(-Q / (R * T_relative))
end

function rate_factor(
    T_relative::T,
    rf::FanHighStrainGSIRateFactor,
) where {T<:Real}
    (; A_0, Q, R) = rf
    return A_0 * exp(-Q / (R * T_relative))
end

function rate_factor(
    T_relative::T,
    arf::ArrheniusRateFactor,
) where {T<:Real}

    (; E_f, T_p1_p2, A_0_p1, A_0_p2, Q_a_p1, Q_a_p2, R) = arf
    if T_relative <= T_p1_p2
        A = E_f * A_0_p1 * exp(-Q_a_p1 / (R * T_relative))
    else
        A = E_f * A_0_p2 * exp(-Q_a_p2 / (R * T_relative))
    end
    return A
end

function rate_factor(
    T_relative::T,
    smr::SmithMorlandRateFactor,
) where {T<:Real}

    (; p1, p2, e1, e2) = smr
    T_bar = (T - T_0) / ΔT
    A = p1 * exp(e1 * T_relative) + p2 * exp(e2 * T_relative)
    return A
end

function rate_factor(
    T_relative::T,
    hrf::HookeRateFactor,
) where {T<:Real}
    (; E_f, A_0, Q_a, C, T_r, k, R) = hrf
    return E_f * A_0 * exp(-Q_a / (R * T_relative) + 3 * C / (T_r - T_relative)^k)
end

function rate_factor(
    T_relative::T,
    ldrf::LliboutryDuvalRateFactor,
) where {T<:Real}
    (; ω, γ, E_f, T_p1_p2, A_0_p1, A_0_p2, Q_a_p1, Q_a_p2, R) = ldrf
    if T_relative <= T_p1_p2
        A_cold = E_f * A_0_p1 * exp(-Q_a_p1 / (R * T_relative))
    else
        A_cold = E_f * A_0_p2 * exp(-Q_a_p2 / (R * T_relative))
    end
    return A_cold * (1 + γ * ω)
end

function rate_factor(
    T_relative::M,
    arf::ARF,
) where {M<:AbstractArray, ARF<:AbstractRateFactor}
    A = similar(T_relative)
    rate_factor!(A, T_relative, arf)
    return A
end

"""
$(TYPEDSIGNATURES)

Get the rate factor `A` based on the temperature relative to the pressure melt point `T_relative` and the rate factor parameterization `arf<:AbstractRateFactor`. Optionally, a mask can be provided to only compute the rate factor for specific indices.
"""
function rate_factor!(
    A,
    T_relative,
    arf::AbstractRateFactor,
)
    map!(x -> rate_factor(x, arf), A, T_relative)
    return
end
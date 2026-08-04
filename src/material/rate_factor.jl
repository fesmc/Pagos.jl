###########################################################
# Structs
###########################################################

"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the rate factor computation via [`rate_factor`](@ref).

!!! note "Rate factors are the one place the time unit lives"
    [`rate_factor`](@ref) returns `A` in `Pa⁻ⁿ yr⁻¹`, matching Pagos's `(m, yr, Pa)`
    base units (see [`Constants`](@ref)). Every other field of every rate factor —
    activation energies, the gas constant, breakpoint temperatures, exponents — is
    time-free, so the pre-exponential factors are the *only* numbers here that change
    if you change the time unit.

    The literature publishes those pre-exponentials in `Pa⁻ⁿ s⁻¹`. Each default below
    is therefore written as `<published value> * SECONDS_PER_YEAR`, so it stays
    diffable against its source. To pass your own values in the published convention,
    give `:second` as the first positional argument:

    ```julia
    ArrheniusRateFactor(:second; A_0_p1 = 3.985e-13, A_0_p2 = 1.916e3)  # Pa⁻³ s⁻¹ in
    ArrheniusRateFactor(A_0_p1 = 1.2578e-5)                             # Pa⁻³ yr⁻¹ in
    ```
"""
abstract type AbstractRateFactor end

"""
$(TYPEDSIGNATURES)

Rate factor for ice viscosity following a constant value.

Glen's flow law with `n = 3`; the default is [cuffey_physics_2010](@citet)'s value for ice
near the melting point, the figure conventionally quoted as `1e-16 Pa⁻³ yr⁻¹`
(`3.2e-24 Pa⁻³ s⁻¹`).

# Fields
- `A::T=1e-16`: rate factor (``\\mathrm{Pa}^{-3}\\,\\mathrm{yr}^{-1}``).
"""
@kwdef struct PrescribedRateFactor{T<:Real} <: AbstractRateFactor
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
 - `A_0_p1`: pre-exponential factor ``A_{0,1}`` (``\\mathrm{Pa}^{-3}\\,\\mathrm{yr}^{-1}``),
   defaulting to Paterson's published ``3.985 \\times 10^{-13}\\,\\mathrm{Pa}^{-3}\\,\\mathrm{s}^{-1}``.
 - `A_0_p2`: pre-exponential factor ``A_{0,2}`` (``\\mathrm{Pa}^{-3}\\,\\mathrm{yr}^{-1}``),
   defaulting to Paterson's published ``1.916 \\times 10^{3}\\,\\mathrm{Pa}^{-3}\\,\\mathrm{s}^{-1}``.
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
@kwdef struct ArrheniusRateFactor{T<:Real} <: AbstractRateFactor
    E_f::T = 1.0
    T_p1_p2::T = 263.15                          # breakpoint temperature following Patterson (1994)
    A_0_p1::T = 3.985e-13 * SECONDS_PER_YEAR     # Pa^-3 s^-1 published -> Pa^-3 yr^-1 internal
    A_0_p2::T = 1.916e3   * SECONDS_PER_YEAR     # piecewise definition following Patterson (1994)
    Q_a_p1::T = 60e3                             # "J mol^-1" activation energy
    Q_a_p2::T = 139e3                            #  piecewise definition following Patterson (1994)
    R::T = 8.314                                 # "J K^-1 mol^-1" universal gas constant
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
@kwdef struct SmithMorlandRateFactor{T<:Real} <: AbstractRateFactor
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
 - `A_0`: pre-exponential factor (``\\mathrm{Pa}^{-3}\\,\\mathrm{yr}^{-1}``), defaulting to
   Hooke's published ``9.302 \\times 10^{-7}\\,\\mathrm{Pa}^{-3}\\,\\mathrm{s}^{-1}``.
 - `Q_a::T=78.8e3`: activation energy (``\\mathrm{J}\\,\\mathrm{mol}^{-1}``).
 - `C::T=0.16612`: proximity-to-melting coefficient (``\\mathrm{K}^{k}``).
 - `T_r::T=273.39`: reference melting temperature (``\\mathrm{K}``).
 - `k::T=1.17`: proximity-to-melting exponent.
 - `R::T=8.314`: universal gas constant (``\\mathrm{J}\\,\\mathrm{K}^{-1}\\,\\mathrm{mol}^{-1}``).
"""
@kwdef struct HookeRateFactor{T<:Real} <: AbstractRateFactor
    E_f::T = 1.0
    A_0::T = 9.302e-7 * SECONDS_PER_YEAR   # Pa^-3 s^-1 published -> Pa^-3 yr^-1 internal
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
 - `A_0_p1`: pre-exponential factor below breakpoint (``\\mathrm{Pa}^{-3}\\,\\mathrm{yr}^{-1}``),
   defaulting to the published ``3.985 \\times 10^{-13}\\,\\mathrm{Pa}^{-3}\\,\\mathrm{s}^{-1}``.
 - `A_0_p2`: pre-exponential factor above breakpoint, same convention.
 - `Q_a_p1::T=60e3`: activation energy below breakpoint (``\\mathrm{J}\\,\\mathrm{mol}^{-1}``).
 - `Q_a_p2::T=139e3`: activation energy above breakpoint.
 - `R::T=8.314`: universal gas constant (``\\mathrm{J}\\,\\mathrm{K}^{-1}\\,\\mathrm{mol}^{-1}``).
"""
@kwdef struct LliboutryDuvalRateFactor{T<:Real} <: AbstractRateFactor
    ω::T = 0.0
    γ::T = 181.25
    E_f::T = 1.0
    T_p1_p2::T = 263.15
    A_0_p1::T = 3.985e-13 * SECONDS_PER_YEAR   # Pa^-3 s^-1 published -> Pa^-3 yr^-1 internal
    A_0_p2::T = 1.916e3   * SECONDS_PER_YEAR
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
 - `A_0`: pre-exponential factor (``\\mathrm{Pa}^{-3.6}\\,\\mathrm{yr}^{-1}``), defaulting to the
   published ``2.951 \\times 10^{-17}\\,\\mathrm{Pa}^{-3.6}\\,\\mathrm{s}^{-1}``.
 - `Q::T=62e3`: activation energy (``\\mathrm{J}\\,\\mathrm{mol}^{-1}``).
 - `R::T=8.314`: universal gas constant (``\\mathrm{J}\\,\\mathrm{K}^{-1}\\,\\mathrm{mol}^{-1}``).
"""
@kwdef struct FanLowStrainGSIRateFactor{T<:Real} <: AbstractRateFactor
    A_0::T = 2.951e-17 * SECONDS_PER_YEAR   # 10^(5.07 - 6*3.6) Pa^-3.6 s^-1 -> Pa^-3.6 yr^-1
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
 - `A_0`: pre-exponential factor (``\\mathrm{Pa}^{-1.9}\\,\\mathrm{m}^{1.2}\\,\\mathrm{yr}^{-1}``),
   defaulting to the published ``4.677 \\times 10^{-13}`` in ``\\mathrm{s}^{-1}``.
 - `Q::T=52e3`: activation energy (``\\mathrm{J}\\,\\mathrm{mol}^{-1}``).
 - `R::T=8.314`: universal gas constant (``\\mathrm{J}\\,\\mathrm{K}^{-1}\\,\\mathrm{mol}^{-1}``).
"""
@kwdef struct FanLowStrainGSS1RateFactor{T<:Real} <: AbstractRateFactor
    A_0::T = 4.677e-13 * SECONDS_PER_YEAR   # 10^(-0.93 - 6*1.9) Pa^-1.9 m^1.2 s^-1 -> yr^-1
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
 - `A_0`: pre-exponential factor (``\\mathrm{Pa}^{-2.5}\\,\\mathrm{m}^{1.9}\\,\\mathrm{yr}^{-1}``),
   defaulting to the published ``4.571 \\times 10^{7}`` in ``\\mathrm{s}^{-1}``.
 - `Q::T=182e3`: activation energy (``\\mathrm{J}\\,\\mathrm{mol}^{-1}``).
 - `R::T=8.314`: universal gas constant (``\\mathrm{J}\\,\\mathrm{K}^{-1}\\,\\mathrm{mol}^{-1}``).
"""
@kwdef struct FanLowStrainGSS2RateFactor{T<:Real} <: AbstractRateFactor
    A_0::T = 4.571e7 * SECONDS_PER_YEAR     # 10^(22.66 - 6*2.5) Pa^-2.5 m^1.9 s^-1 -> yr^-1
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
 - `A_0`: pre-exponential factor (``\\mathrm{Pa}^{-3.5}\\,\\mathrm{yr}^{-1}``), defaulting to the
   published ``7.943 \\times 10^{-10}\\,\\mathrm{Pa}^{-3.5}\\,\\mathrm{s}^{-1}``.
 - `Q::T=90e3`: activation energy (``\\mathrm{J}\\,\\mathrm{mol}^{-1}``).
 - `R::T=8.314`: universal gas constant (``\\mathrm{J}\\,\\mathrm{K}^{-1}\\,\\mathrm{mol}^{-1}``).
"""
@kwdef struct FanHighStrainGSIRateFactor{T<:Real} <: AbstractRateFactor
    A_0::T = 7.943e-10 * SECONDS_PER_YEAR   # 10^(11.9 - 6*3.5) Pa^-3.5 s^-1 -> Pa^-3.5 yr^-1
    Q::T = 90e3
    R::T = 8.314
end

###########################################################
# Time-unit handling
###########################################################

# Which fields carry the time unit, i.e. which ones `time_unit` rescales. Everything not
# listed is time-free and passes through untouched. `SmithMorlandRateFactor` falls back to
# the empty default because its coefficients are dimensionless — that law's time unit sits
# in `SmithMorlandCreep`'s `D_0` instead.
_prefactor_fields(::Type{<:AbstractRateFactor})         = ()
_prefactor_fields(::Type{<:PrescribedRateFactor})       = (:A,)
_prefactor_fields(::Type{<:ArrheniusRateFactor})        = (:A_0_p1, :A_0_p2)
_prefactor_fields(::Type{<:LliboutryDuvalRateFactor})   = (:A_0_p1, :A_0_p2)
_prefactor_fields(::Type{<:HookeRateFactor})            = (:A_0,)
_prefactor_fields(::Type{<:FanLowStrainGSIRateFactor})  = (:A_0,)
_prefactor_fields(::Type{<:FanLowStrainGSS1RateFactor}) = (:A_0,)
_prefactor_fields(::Type{<:FanLowStrainGSS2RateFactor}) = (:A_0,)
_prefactor_fields(::Type{<:FanHighStrainGSIRateFactor}) = (:A_0,)

_time_unit_scale(u::Symbol) =
    u === :year   ? 1.0 :
    u === :second ? SECONDS_PER_YEAR :
    throw(ArgumentError("time_unit must be :year or :second, got :$u"))

"""
$(TYPEDSIGNATURES)

Construct a rate factor from pre-exponential factors given in `time_unit` (`:second` for
values as published in SI, `:year` for Pagos's internal convention — see
[`AbstractRateFactor`](@ref)).

Only the pre-exponential fields are rescaled; activation energies, temperatures and the
gas constant are time-free and pass through. Keywords you do not pass keep their defaults,
which are already internal, so partial overrides are safe:

```julia
ArrheniusRateFactor(:second; A_0_p1 = 3.985e-13)   # A_0_p2 keeps its default
```
"""
function (::Type{RF})(time_unit::Symbol; kwargs...) where {RF<:AbstractRateFactor}
    return RF(; _rescaled_kwargs(RF, time_unit, kwargs)...)
end

# `PrescribedRateFactor` is the only rate factor with exactly one field, so `@kwdef` also
# generates an *untyped* positional constructor `PrescribedRateFactor{T}(A)`. That matches
# a lone `Symbol` whatever the bound on `T`, making it ambiguous with the generic method
# above; this more specific method resolves it. Aqua's ambiguity check guards the case.
function (::Type{PrescribedRateFactor{T}})(time_unit::Symbol; kwargs...) where {T<:Real}
    return PrescribedRateFactor{T}(; _rescaled_kwargs(PrescribedRateFactor, time_unit, kwargs)...)
end

function _rescaled_kwargs(::Type{RF}, time_unit::Symbol, kwargs) where {RF<:AbstractRateFactor}
    scale = _time_unit_scale(time_unit)
    fields = _prefactor_fields(RF)
    return NamedTuple(k => (k in fields ? v * scale : v) for (k, v) in pairs(kwargs))
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
    crf::PrescribedRateFactor,
) where {T<:Real}
    return crf.A
end

function rate_factor(
    T_relative::T,
    rf::F,
) where {T<:Real, F<:AbstractRateFactor}
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
function rate_factor!(A::AbstractArray, T_relative, arf::AbstractRateFactor)
    pointwise!(rate_factor, A, (T_relative,), (arf,))
    return nothing
end
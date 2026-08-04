###########################################################
# Structs
###########################################################

"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the creep function computation via [`creep`](@ref) and [`creep!`](@ref).
"""
abstract type AbstractCreep end

"""
$(TYPEDSIGNATURES)

Regularized Glen-Nye creep function following [glen_creep_1955](@citet) and [nye_distribution_1957](@citet):

```math
\\begin{aligned}
f(\\sigma_e) = \\sigma_e^{n - 1} + \\sigma_0^{n - 1}
\\end{aligned}
```

The floor ``\\sigma_0^{n-1}`` prevents the viscosity singularity at zero stress.
In practice ``\\sigma_0 \\ll \\sigma_e`` for realistic stress levels, so the
regularization is only active near ice divides where stresses approach zero.

# Fields
 - `n::T=3.0`: Glen-Nye's flow law exponent.
 - `σ_0::T=1e-6`: regularization stress (``\\mathrm{Pa}``).
"""
@kwdef struct GlenNyeCreep{T} <: AbstractCreep
    n::T = 3.0
    σ_0::T = 1e-6
end

"""
$(TYPEDSIGNATURES)

Polynomial creep function following [smith_viscous_1981](@citet):

```math
\\begin{aligned}
f(\\sigma_e) = \\frac{D_0}{\\sigma_0}
\\left( p_0 + p_2 \\left(\\frac{\\sigma_e}{\\sigma_0}\\right)^{\\!2}
           + p_4 \\left(\\frac{\\sigma_e}{\\sigma_0}\\right)^{\\!4} \\right)
\\end{aligned}
```

The polynomial is fit to laboratory data. Unlike the Glen-Nye family, this
parameterization is calibrated to pair with [`SmithMorlandRateFactor`](@ref)
via [`SmithMorlandFlowLaw`](@ref); see that docstring for units and caveats.

!!! warning
    The field names `p1`, `p2`, `p3` in the docstring differ from the actual struct
    fields `p0`, `p2`, `p4` — use the struct fields directly.

# Fields
 - `p0::T=0.3336`: zeroth-order polynomial coefficient.
 - `p2::T=0.3200`: second-order polynomial coefficient.
 - `p4::T=0.02963`: fourth-order polynomial coefficient.
 - `D_0`: strain-rate pre-factor (``\\mathrm{Pa}^{-1}\\,\\mathrm{yr}^{-1}``), defaulting to the
   published ``3.169 \\times 10^{-8}\\,\\mathrm{Pa}^{-1}\\,\\mathrm{s}^{-1}``. That figure is
   ``1/(\\mathrm{s}\\,\\mathrm{yr}^{-1})`` to four digits, i.e. the original calibration was
   itself per-year and had been divided by it — so the internal value comes back to `1.0`.
 - `σ_0::T=1e5`: reference stress (``\\mathrm{Pa}``).
"""
@kwdef struct SmithMorlandCreep{T} <: AbstractCreep
    p0::T = 0.3336
    p2::T = 0.3200
    p4::T = 0.02963
    D_0::T = 3.169e-8 * SECONDS_PER_YEAR   # Pa^-1 s^-1 published -> Pa^-1 yr^-1 (= 1.0)
    σ_0::T = 1e5
end

"""
$(TYPEDSIGNATURES)

Composite creep function following [pettit_ice_2003](@citet), combining a linear
(Newtonian) and a nonlinear (Glen-type) contribution:

```math
\\begin{aligned}
f(\\sigma_e) = A_{\\mathrm{lin}} + A_{\\mathrm{nl}} \\, \\sigma_e^{n - 1}
\\end{aligned}
```

At low stresses the linear term dominates and the viscosity approaches the
Newtonian limit; at high stresses the Glen power-law term takes over.
Reduces to [`GlenNyeCreep`](@ref) when ``A_{\\mathrm{lin}} = 0``.
Suitable for modelling flow near ice divides.

# Fields
 - `n::T=3.0`: Glen-type flow exponent.
 - `A_lin::T=0.0`: linear (Newtonian) coefficient (``\\mathrm{Pa}^{n-1}``).
 - `A_nl::T=1.0`: nonlinear (Glen) coefficient (dimensionless).
"""
@kwdef struct PettitWaddingtonCreep{T} <: AbstractCreep
    n::T = 3.0
    A_lin::T = 0.0
    A_nl::T = 1.0
end

"""
$(TYPEDSIGNATURES)

Composite creep function following [goldsby_superplastic_2001](@citet). Diffusion
creep, grain-boundary sliding (GBS), and basal dislocation slip operate in parallel,
with GBS and basal slip coupled in series:

```math
\\begin{aligned}
\\dot{\\varepsilon}_{\\mathrm{diff}}  &= A_{\\mathrm{diff}} \\, d^{-2} \\, \\sigma_e \\\\
\\dot{\\varepsilon}_{\\mathrm{gbs}}   &= A_{\\mathrm{gbs}}  \\, d^{-p_{\\mathrm{gbs}}} \\, \\sigma_e^{n_{\\mathrm{gbs}}} \\\\
\\dot{\\varepsilon}_{\\mathrm{basal}} &= A_{\\mathrm{basal}} \\, \\sigma_e^{n_{\\mathrm{basal}}} \\\\[4pt]
\\dot{\\varepsilon}_{\\mathrm{tot}}   &= \\dot{\\varepsilon}_{\\mathrm{diff}}
    + \\left(\\dot{\\varepsilon}_{\\mathrm{gbs}}^{-1} + \\dot{\\varepsilon}_{\\mathrm{basal}}^{-1}\\right)^{-1}
\\end{aligned}
```

The creep function returned is ``f(\\sigma_e) = \\dot{\\varepsilon}_{\\mathrm{tot}} / \\sigma_e``.
The pre-factors ``A_{\\mathrm{diff}}``, ``A_{\\mathrm{gbs}}``, ``A_{\\mathrm{basal}}``
carry their own Arrhenius temperature dependence and must be evaluated at the local
temperature before constructing this struct. They are caller-supplied, so the caller also
owns the conversion: pass them in ``\\mathrm{yr}^{-1}``, not the ``\\mathrm{s}^{-1}`` of
the source literature (see [`AbstractRateFactor`](@ref)).

# Fields
 - `d::T`: mean grain size (``\\mathrm{m}``).
 - `A_diff::T`: diffusion creep pre-factor (``\\mathrm{Pa}^{-1}\\,\\mathrm{m}^2\\,\\mathrm{yr}^{-1}``).
 - `A_gbs::T`: GBS pre-factor (``\\mathrm{Pa}^{-n_{\\mathrm{gbs}}}\\,\\mathrm{m}^{p_{\\mathrm{gbs}}}\\,\\mathrm{yr}^{-1}``).
 - `A_basal::T`: basal slip pre-factor (``\\mathrm{Pa}^{-n_{\\mathrm{basal}}}\\,\\mathrm{yr}^{-1}``).
 - `p_gbs::T=1.4`: grain-size exponent for GBS.
 - `n_gbs::T=1.8`: stress exponent for GBS.
 - `n_basal::T=2.4`: stress exponent for basal slip.
"""
@kwdef struct GoldsbyKohlstedtCreep{T} <: AbstractCreep
    d::T
    A_diff::T
    A_gbs::T
    A_basal::T
    p_gbs::T = 1.4
    n_gbs::T = 1.8
    n_basal::T = 2.4
end

"""
$(TYPEDSIGNATURES)

The full four-component Goldsby–Kohlstedt creep function of
[goldsby_superplastic_2001](@citet), matching PISM's `gk` flow law. Diffusional flow,
dislocation creep, grain-boundary sliding (GBS) and basal (easy) slip, with GBS and basal
slip coupled in series and the other two in parallel:

```math
\\begin{aligned}
\\dot{\\varepsilon}_{\\mathrm{diff}}  &= \\frac{42\\, V_m}{R T d^{2}}
    \\left( D_v + \\frac{\\pi \\delta D_b}{d} \\right) \\sigma_e \\\\
\\dot{\\varepsilon}_{\\mathrm{disl}}  &= A_{\\mathrm{disl}} \\, \\sigma_e^{n_{\\mathrm{disl}}}
    \\, e^{-(Q_{\\mathrm{disl}} + pV)/RT} \\\\
\\dot{\\varepsilon}_{\\mathrm{basal}} &= A_{\\mathrm{basal}} \\, \\sigma_e^{n_{\\mathrm{basal}}}
    \\, e^{-(Q_{\\mathrm{basal}} + pV)/RT} \\\\
\\dot{\\varepsilon}_{\\mathrm{gbs}}   &= A_{\\mathrm{gbs}} \\, \\sigma_e^{n_{\\mathrm{gbs}}}
    \\, d^{-p_{\\mathrm{gbs}}} \\, e^{-(Q_{\\mathrm{gbs}} + pV)/RT} \\\\[4pt]
\\dot{\\varepsilon}_{\\mathrm{tot}}   &= \\dot{\\varepsilon}_{\\mathrm{diff}}
    + \\dot{\\varepsilon}_{\\mathrm{disl}}
    + \\left(\\dot{\\varepsilon}_{\\mathrm{basal}}^{-1}
          + \\dot{\\varepsilon}_{\\mathrm{gbs}}^{-1}\\right)^{-1}
\\end{aligned}
```

with ``D_v = D_{0v} e^{-Q_v/RT}`` and ``D_b = D_{0b} e^{-Q_b/RT}``. The creep function
returned is ``f(\\sigma_e) = \\dot{\\varepsilon}_{\\mathrm{tot}} / \\sigma_e``.

!!! warning "Untested"
    Nothing in `test/` exercises this yet — it is transcribed from PISM's
    `src/rheology/GoldsbyKohlstedt.cc` and checked only by eye. Validate before relying on
    it quantitatively.

Three things distinguish it from the simplified [`GoldsbyKohlstedtCreep`](@ref): the extra
dislocation-creep component (``n = 4``), cold/warm Arrhenius branches for dislocation creep
and GBS, and a pressure dependence entering through the activation volume ``V`` as
``e^{-(Q + pV)/RT}``. Because of that pressure and temperature dependence this struct is
evaluated at a point rather than taking pre-evaluated pre-factors — pass the local
`temperature` and `pressure` and construct it pointwise.

!!! note "`temperature` is pressure-adjusted"
    Give the temperature relative to the pressure melting point (see
    [`LinearPressureMeltingPoint`](@ref)), the same convention every
    [`AbstractRateFactor`](@ref) here uses — *not* the in-situ temperature. `pressure` is
    separate and feeds only the ``pV`` activation-volume term.

# Fields
 - `temperature::T`: pressure-adjusted temperature (``\\mathrm{K}``).
 - `pressure::T`: ice overburden pressure (``\\mathrm{Pa}``).
 - `d::T=1.0e-3`: mean grain size (``\\mathrm{m}``).
 - `R::T=8.314`: universal gas constant (``\\mathrm{J}\\,\\mathrm{K}^{-1}\\,\\mathrm{mol}^{-1}``).
 - `V_act::T=-13.0e-6`: activation volume (``\\mathrm{m}^3\\,\\mathrm{mol}^{-1}``). Negative,
   so pressure *softens* the ice.
 - `disl_T_crit::T=258.0`, `disl_A_cold`, `disl_A_warm`, `disl_n::T=4.0`,
   `disl_Q_cold::T=60e3`, `disl_Q_warm::T=180e3`: dislocation creep. Pre-factors are
   ``\\mathrm{Pa}^{-4}\\,\\mathrm{yr}^{-1}``, from the published
   ``4.0\\times10^{-19}`` / ``6.0\\times10^{4}\\,\\mathrm{Pa}^{-4}\\,\\mathrm{s}^{-1}``.
 - `gbs_T_crit::T=255.0`, `gbs_A_cold`, `gbs_A_warm`, `gbs_n::T=1.8`,
   `gbs_Q_cold::T=49e3`, `gbs_Q_warm::T=192e3`, `gbs_p::T=1.4`: grain-boundary sliding.
 - `basal_A`, `basal_n::T=2.4`, `basal_Q::T=60e3`: basal (easy) slip.
 - `diff_T_crit::T=258.0`, `diff_V_m::T=1.97e-5`, `diff_D_0v`, `diff_Q_v::T=59.4e3`,
   `diff_D_0b`, `diff_Q_b::T=49e3`, `diff_delta::T=9.04e-10`: diffusional flow. Above
   `diff_T_crit` the grain-boundary diffusivity is multiplied by `diff_coble_factor`.
 - `diff_coble_factor::T=1000.0`: Coble-creep scaling applied to ``D_b`` when warm.
"""
@kwdef struct GoldsbyKohlstedt4Creep{T<:Real} <: AbstractCreep
    temperature::T
    pressure::T
    d::T = 1.0e-3
    R::T = 8.314
    V_act::T = -13.0e-6

    ## dislocation creep
    disl_T_crit::T = 258.0
    disl_A_cold::T = 4.0e-19 * SECONDS_PER_YEAR   # Pa^-4 s^-1 published -> Pa^-4 yr^-1
    disl_A_warm::T = 6.0e4 * SECONDS_PER_YEAR
    disl_n::T = 4.0
    disl_Q_cold::T = 60.0e3
    disl_Q_warm::T = 180.0e3

    ## grain-boundary sliding
    gbs_T_crit::T = 255.0
    gbs_A_cold::T = 6.1811e-14 * SECONDS_PER_YEAR # Pa^-1.8 m^1.4 s^-1 -> ... yr^-1
    gbs_A_warm::T = 4.7547e15 * SECONDS_PER_YEAR
    gbs_n::T = 1.8
    gbs_Q_cold::T = 49.0e3
    gbs_Q_warm::T = 192.0e3
    gbs_p::T = 1.4

    ## basal (easy) slip
    basal_A::T = 2.1896e-7 * SECONDS_PER_YEAR     # Pa^-2.4 s^-1 -> Pa^-2.4 yr^-1
    basal_n::T = 2.4
    basal_Q::T = 60.0e3

    ## diffusional flow
    diff_T_crit::T = 258.0
    diff_V_m::T = 1.97e-5
    diff_D_0v::T = 9.10e-4 * SECONDS_PER_YEAR     # m^2 s^-1 -> m^2 yr^-1
    diff_Q_v::T = 59.4e3
    diff_D_0b::T = 5.8e-4 * SECONDS_PER_YEAR
    diff_Q_b::T = 49.0e3
    diff_delta::T = 9.04e-10
    diff_coble_factor::T = 1000.0
end

"""
$(TYPEDSIGNATURES)

Three-component low-strain creep function following [fan_flow_2025](@citet),
summing contributions from one grain-size insensitive (GSI) dislocation-creep
component and two grain-size sensitive (GSS) disGBS components:

```math
\\begin{aligned}
\\dot{\\varepsilon}_{\\mathrm{GSI}}  &= A_{\\mathrm{GSI}}  \\, \\sigma_e^{n_{\\mathrm{GSI}}} \\\\
\\dot{\\varepsilon}_{\\mathrm{GSS1}} &= A_{\\mathrm{GSS1}} \\, \\sigma_e^{n_{\\mathrm{GSS1}}} \\, d^{-p_{\\mathrm{GSS1}}} \\\\
\\dot{\\varepsilon}_{\\mathrm{GSS2}} &= A_{\\mathrm{GSS2}} \\, \\sigma_e^{n_{\\mathrm{GSS2}}} \\, d^{-p_{\\mathrm{GSS2}}} \\\\[4pt]
f(\\sigma_e) &= \\frac{\\dot{\\varepsilon}_{\\mathrm{GSI}} + \\dot{\\varepsilon}_{\\mathrm{GSS1}} + \\dot{\\varepsilon}_{\\mathrm{GSS2}}}{\\sigma_e}
\\end{aligned}
```

The two GSS components represent a single disGBS mechanism whose apparent
activation energy increases near the pressure melting point, capturing the
continuous temperature dependence of ice viscosity without a threshold
discontinuity.

The pre-factors ``A_{\\mathrm{GSI}}``, ``A_{\\mathrm{GSS1}}``, ``A_{\\mathrm{GSS2}}``
already encode the Arrhenius temperature dependence and must be evaluated at the
local temperature via [`FanLowStrainGSIRateFactor`](@ref),
[`FanLowStrainGSS1RateFactor`](@ref), and [`FanLowStrainGSS2RateFactor`](@ref)
before constructing this struct. Use [`FanLowStrainFlowLaw`](@ref) for a
convenience constructor that handles this step automatically.

# Fields
 - `d::T`: mean grain size (``\\mathrm{m}``).
 - `A_GSI::T`: temperature-evaluated GSI pre-factor (``\\mathrm{Pa}^{-n_{\\mathrm{GSI}}}\\,\\mathrm{yr}^{-1}``).
 - `A_GSS1::T`: temperature-evaluated GSS1 pre-factor (``\\mathrm{Pa}^{-n_{\\mathrm{GSS1}}}\\,\\mathrm{m}^{p_{\\mathrm{GSS1}}}\\,\\mathrm{yr}^{-1}``).
 - `A_GSS2::T`: temperature-evaluated GSS2 pre-factor (``\\mathrm{Pa}^{-n_{\\mathrm{GSS2}}}\\,\\mathrm{m}^{p_{\\mathrm{GSS2}}}\\,\\mathrm{yr}^{-1}``).
 - `n_GSI::T=3.6`: stress exponent for the GSI component.
 - `n_GSS1::T=1.9`: stress exponent for the GSS1 component.
 - `n_GSS2::T=2.5`: stress exponent for the GSS2 component.
 - `p_GSS1::T=1.2`: grain-size exponent for the GSS1 component.
 - `p_GSS2::T=1.9`: grain-size exponent for the GSS2 component.
"""
@kwdef struct FanLowStrainCreep{T} <: AbstractCreep
    d::T
    A_GSI::T
    A_GSS1::T
    A_GSS2::T
    n_GSI::T = 3.6
    n_GSS1::T = 1.9
    n_GSS2::T = 2.5
    p_GSS1::T = 1.2
    p_GSS2::T = 1.9
end

###########################################################
# Dispatch
###########################################################

"""
$(TYPEDSIGNATURES)

Compute the creep function value based on the effective stress `σ_e` and the creep function parameterization `law<:AbstractCreep`.
"""
function creep(
    σ_e::T,
    law::SmithMorlandCreep,
) where {T<:Real}
    (; p0, p2, p4, D_0, σ_0) = law
    return D_0 / σ_0 * (p0 + p2 * (σ_e / σ_0)^2 + p4 * (σ_e / σ_0)^4)
end

function creep(
    σ_e::T,
    law::GlenNyeCreep,
) where {T<:Real}
    (; n, σ_0) = law
    return σ_e^(n - 1) + σ_0^(n - 1)
end

function creep(
    σ_e::T,
    law::PettitWaddingtonCreep,
) where {T<:Real}
    (; n, A_lin, A_nl) = law
    return A_lin + A_nl * σ_e^(n - 1)
end

function creep(
    σ_e::T,
    law::GoldsbyKohlstedtCreep,
) where {T<:Real}
    (; d, A_diff, A_gbs, A_basal, p_gbs, n_gbs, n_basal) = law
    ε̇_diff  = A_diff * d^(-2) * σ_e
    ε̇_gbs   = A_gbs * d^(-p_gbs) * σ_e^n_gbs
    ε̇_basal = A_basal * σ_e^n_basal
    ε̇_eff   = ε̇_diff + inv(inv(ε̇_gbs) + inv(ε̇_basal))
    return ε̇_eff / σ_e
end

function creep(
    σ_e::T,
    law::GoldsbyKohlstedt4Creep,
) where {T<:Real}
    (; temperature, pressure, d, R, V_act) = law
    RT = R * temperature
    pV = pressure * V_act

    ## Diffusional flow. Newtonian, so ε̇ ∝ σ_e and this contributes a stress-independent
    ## term to f. Coble creep enhances grain-boundary diffusion above `diff_T_crit`.
    D_v = law.diff_D_0v * exp(-law.diff_Q_v / RT)
    D_b = law.diff_D_0b * exp(-law.diff_Q_b / RT)
    D_b = temperature > law.diff_T_crit ? D_b * law.diff_coble_factor : D_b
    f_diff = 42 * law.diff_V_m * (D_v + π * law.diff_delta * D_b / d) / (RT * d^2)

    A_disl, Q_disl = temperature > law.disl_T_crit ?
        (law.disl_A_warm, law.disl_Q_warm) : (law.disl_A_cold, law.disl_Q_cold)
    f_disl = A_disl * σ_e^(law.disl_n - 1) * exp(-(Q_disl + pV) / RT)

    f_basal = law.basal_A * σ_e^(law.basal_n - 1) * exp(-(law.basal_Q + pV) / RT)

    A_gbs, Q_gbs = temperature > law.gbs_T_crit ?
        (law.gbs_A_warm, law.gbs_Q_warm) : (law.gbs_A_cold, law.gbs_Q_cold)
    f_gbs = A_gbs * σ_e^(law.gbs_n - 1) * d^(-law.gbs_p) * exp(-(Q_gbs + pV) / RT)

    # Series coupling as `inv(inv + inv)` rather than the algebraically equal
    # `f_basal*f_gbs/(f_basal+f_gbs)`: at σ_e = 0 both terms vanish and the product form
    # gives 0/0, while this one gives inv(Inf) = 0. Same idiom as `GoldsbyKohlstedtCreep`.
    return f_diff + f_disl + inv(inv(f_basal) + inv(f_gbs))
end

function creep(
    σ_e::T,
    law::FanLowStrainCreep,
) where {T<:Real}
    (; d, A_GSI, A_GSS1, A_GSS2, n_GSI, n_GSS1, n_GSS2, p_GSS1, p_GSS2) = law
    ε̇_GSI  = A_GSI  * σ_e^n_GSI
    ε̇_GSS1 = A_GSS1 * σ_e^n_GSS1 * d^(-p_GSS1)
    ε̇_GSS2 = A_GSS2 * σ_e^n_GSS2 * d^(-p_GSS2)
    return (ε̇_GSI + ε̇_GSS1 + ε̇_GSS2) / σ_e
end

function creep(
    σ_e::M,
    law,
) where {M<:AbstractArray}
    cf = similar(σ_e)
    creep!(cf, σ_e, law)
    return cf
end

"""
$(TYPEDSIGNATURES)

Same as [`creep`](@ref) but operates in place.
"""
function creep!(cf::AbstractArray, σ_e, law)
    pointwise!(creep, cf, (σ_e,), (law,))
    return nothing
end
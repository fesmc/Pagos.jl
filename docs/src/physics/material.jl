#=

# [Material](@id material)

The deformation of polycrystalline ice under shear stress ``\tau`` can be described by a flow law that depends on the ice viscosity ``\eta`` and yields the shear (deformation) rate ``\dot{\gamma}``, which can be expressed as:

```math
\begin{aligned}
\dot{\gamma} = \dfrac{1}{\eta(T, p, \tau)} \, \tau
\end{aligned}
```

The viscosity itself depends on temperature (``T``), pressure (``p``), and shear stress (``\tau``). It can be factorized as:

```math
\begin{aligned}
\eta(T, p, \tau) = \frac{1}{2 \, A(T, p) f(|\tau|)}
\end{aligned}
```

where ``A(T, p)`` is the rate factor, and ``f(|\tau|)`` is the creep function. Let's first look at the rate factor.

## [Rate Factor](@id rate_factor)

The rate factor ``A(T, p)`` captures the dependence of ice deformation on temperature and pressure. It typically increases with temperature and decreases with pressure, reflecting the fact that warmer ice deforms more easily, while higher pressure tends to inhibit deformation.

The rate factor ``A(T, p)`` is usually modeled by an [`ArrheniusRateFactor`](@ref), but any other subtype of [`AbstractRateFactor`](@ref) can be used. Here we stick to this standard choice to illustrate the typical temperature dependence of the rate factor:

=#

using Pagos, CairoMakie
set_theme!(theme_latexfonts())
arrhenius_rate_factor = ArrheniusRateFactor()
T_relative_kelvin = range(-50, stop = 10, step = 0.1) .+ 273.15
A = rate_factor(T_relative_kelvin, arrhenius_rate_factor)
fig = plot_rate_factor(T_relative_kelvin .- 273.15, A)

#=
This matches Fig. 4.5 of the holy bible, a.k.a [greve_dynamics_2009](@citet)!

## [Make your own rate factor model](@id custom_rate_factor)

Besides the many choices provided by the subtypes of [`AbstractRateFactor`](@ref), you can also implement your own rate factor models. For instance, one could smooth the kink in the Arrhenius law at the breakpoint temperature by using a linear combination of the two Arrhenius terms that relies on a sigmoid function:
=#

@kwdef struct SmoothArrheniusRateFactor{T} <: AbstractRateFactor
    arr::ArrheniusRateFactor{T} = ArrheniusRateFactor{T}()  # underlying Arrhenius parameters
    β::T = 1.0                                              # smoothing exponent
end

function Pagos.rate_factor(T_relative::T, rf::SmoothArrheniusRateFactor) where {T<:Real}
    (; E_f, T_p1_p2, A_0_p1, A_0_p2, Q_a_p1, Q_a_p2, R) = rf.arr
    α = 1 / (1 + exp(-rf.β * (T_relative - T_p1_p2)))   # smooth transition function
    A = E_f * (A_0_p1 * exp(-Q_a_p1 / (R * T_relative))) ^ (1 - α) *
        (A_0_p2 * exp(-Q_a_p2 / (R * T_relative))) ^ α
    return A
end

smooth_rf = SmoothArrheniusRateFactor{Float64}()
A = rate_factor(T_relative_kelvin, smooth_rf)
fig_smooth = plot_rate_factor(T_relative_kelvin .- 273.15, A)

#=

If you got lost in mathematical details, worry not! The key point is that you can implement any rate factor model you like by defining a new subtype of [`AbstractRateFactor`](@ref) and implementing the `rate_factor` method for it. The same applies to any abstract type defined in the [`API reference`](@ref) and allows the user to create custom models.

## [Creep function](@id creep_function)

The creep function ``f(|\tau|)`` describes how the viscosity changes with applied shear stress. A commonly used creep function is [`GlenNyeCreep`](@ref) [glen_creep_1955, nye_distribution_1957](@citep), which we here illustrate with a regularized version that prevents a viscosity singularity at zero stress and stress values that typically arise in ice sheets:
=#

glennye_creep = GlenNyeCreep()
σ_e = range(0, stop = 100, step = 0.1) .* 1f3
cf = creep(σ_e, glennye_creep)
fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = L"Effective stress $\sigma_e$ (kPa)",
    ylabel = L"Creep function $f(\sigma_e)$ (Pa²)",
)
lines!(ax, σ_e ./ 1f3, cf)
fig

#=

## Comparing all parameterizations

When choosing between different parameterizations, users might struggle in understanding the differences and implications of each choice. To help with this, we provide a comparison figure that places all built-in rate factor and creep function parameterizations side by side.

**Left panel — rate factor ``A(T)``** (log scale):
[`ArrheniusRateFactor`](@ref) implements the piecewise Arrhenius law of
[paterson_physics_1994](@citet) and exhibits a kink at the breakpoint temperature
``T^{*} = -10°\mathrm{C}`` where the two sets of activation parameters switch.
[`HookeRateFactor`](@ref) is a continuous alternative that augments the standard
Arrhenius expression with a proximity-to-melting correction: the term
``3C\,(T_r - T)^{-k}`` produces a steep upturn as temperature approaches the pressure
melting point, capturing the observed rapid softening of near-temperate ice without an
artificial kink.
[`LliboutryDuvalRateFactor`](@ref) extends the piecewise Arrhenius law with a
water-content enhancement ``(1 + \gamma\,\omega)``: at ``\omega = 0`` it is identical to
[`ArrheniusRateFactor`](@ref), and increasing ``\omega`` shifts the entire curve upward
proportionally.

**Right panel — creep function ``f(\sigma_e)``** (log–log scale):
[`GlenNyeCreep`](@ref) with and without regularization ``\sigma_0^{n-1}`` that prevents a
viscosity singularity at zero stress; here ``\sigma_0 = 1\,\mathrm{kPa}`` is chosen to
make the transition visible, whereas the default ``\sigma_0 = 10^{-6}\,\mathrm{Pa}`` is
negligible at realistic stress levels.
[`PettitWaddingtonCreep`](@ref) adds a linear (Newtonian) term ``A_\mathrm{lin}`` that
dominates at low stresses near ice divides, recovering Glen-type power-law behavior at
high stresses; the two curves illustrate how the transition stress ``\sigma_t \approx
\sqrt{A_\mathrm{lin}}`` shifts with the choice of ``A_\mathrm{lin}``.

!!! note "Smith–Morland and Goldsby–Kohlstedt parameterizations"
    [`SmithMorlandCreep`](@ref) and [`SmithMorlandRateFactor`](@ref) use a different
    convention: the rate factor is dimensionless (input is a normalized temperature
    ``\bar{T} = (T - T_0)/\Delta T``) and the creep function carries units
    ``\mathrm{Pa^{-1}\,s^{-1}}``, unlike the ``\mathrm{Pa^{n-1}}`` units of the Glen-Nye
    family. They are designed to pair together via [`SmithMorlandFlowLaw`](@ref) and
    cannot be compared directly on the axes above.
    [`GoldsbyKohlstedtCreep`](@ref) requires temperature-dependent pre-factors
    ``A_\mathrm{diff}``, ``A_\mathrm{gbs}``, ``A_\mathrm{basal}`` that carry their own
    Arrhenius temperature dependence and must be evaluated at the temperature of interest
    before being passed to [`creep`](@ref).

=#

T_K = range(223.15, stop = 273.15, length = 300)
T_C = T_K .- 273.15
A_arrhenius = rate_factor(T_K, ArrheniusRateFactor())
A_hooke     = rate_factor(T_K, HookeRateFactor())
A_ld1       = rate_factor(T_K, LliboutryDuvalRateFactor(ω = 0.01))
A_ld5       = rate_factor(T_K, LliboutryDuvalRateFactor(ω = 0.05))

σ_range = 10 .^ range(2.0, log10(5e5), length = 300)   # 0.1 kPa – 500 kPa
cf_glen    = creep(σ_range, GlenNyeCreep(σ_0 = 0.0))   # unregularized
cf_reg     = creep(σ_range, GlenNyeCreep(σ_0 = 1e3))
cf_pw_low  = creep(σ_range, PettitWaddingtonCreep(A_lin = 1e6))
cf_pw_high = creep(σ_range, PettitWaddingtonCreep(A_lin = 1e8))

fig_cmp = Figure(size = (900, 420))

ax_rf = Axis(fig_cmp[1, 1],
    xlabel = "Temperature (°C)",
    ylabel = L"Rate factor $A$ (Pa$^{-3}$ s$^{-1}$)",
    yscale = log10,
    title  = "Rate factor parameterizations",
)
lines!(ax_rf, T_C, A_arrhenius, label = "Arrhenius (Paterson 1994)")
lines!(ax_rf, T_C, A_hooke,     label = "Hooke (1981)")
lines!(ax_rf, T_C, A_ld1,       label = L"Lliboutry-Duval ($\omega = 0.01$)")
lines!(ax_rf, T_C, A_ld5,       label = L"Lliboutry-Duval ($\omega = 0.05$)")
axislegend(ax_rf, position = :lt, labelsize = 11)

ax_cf = Axis(fig_cmp[1, 2],
    xlabel = L"Effective stress $\sigma_e$ (kPa)",
    ylabel = L"Creep function $f(\sigma_e)$ (Pa$^2$)",
    yscale = log10,
    xscale = log10,
    title  = "Creep function parameterizations",
)
lines!(ax_cf, σ_range ./ 1e3, cf_glen,    label = L"Glen-Nye ($n = 3$)", linestyle = :dash)
lines!(ax_cf, σ_range ./ 1e3, cf_reg,     label = L"Regularized Glen-Nye ($\sigma_0 = 1$ kPa)")
lines!(ax_cf, σ_range ./ 1e3, cf_pw_low,  label = L"Pettit-Waddington ($A_\mathrm{lin} = 10^6$ Pa$^2$)", linestyle = :dash)
lines!(ax_cf, σ_range ./ 1e3, cf_pw_high, label = L"Pettit-Waddington ($A_\mathrm{lin} = 10^8$ Pa$^2$)")
axislegend(ax_cf, position = :rb, labelsize = 11)

fig_cmp

#=
## [Flow law](@id flow_law)

By combining the rate factor and creep function, we can compute the ice viscosity using the [`RateCreepFlowLaw`](@ref):
=#

flowlaw = RateCreepFlowLaw(arrhenius_rate_factor, glennye_creep)
T = [0, -10, -20]
A = rate_factor(T .+ 273.15, arrhenius_rate_factor)
η = [viscosity(a, cf, flowlaw) for a in A]
fig = plot_ice_viscosity(η, σ_e, T)

#=

### Fan et al. (2025) multicomponent flow laws

[fan_flow_2025](@citet) constrain new flow laws from Bayesian inference on 70 years
of laboratory data. Two are provided:

- **[`FanLowStrainFlowLaw`](@ref)** — three-component law (grain-size insensitive
  dislocation creep + two grain-size sensitive disGBS components) calibrated on
  low-strain (1–2 %) data. It requires a mean grain size ``d`` and a temperature ``T``.
- **[`FanHighStrainFlowLaw`](@ref)** — single GSI component with ``n = 3.5`` and
  ``Q = 90`` kJ mol⁻¹, calibrated on tertiary-creep / flow-stress data. Suitable for
  large-scale ice-sheet models where ice has reached microstructural steady state.

**Left panel** — viscosity versus effective stress at ``T = -5`` °C. The low-strain
law is shown for two grain sizes; its effective stress exponent varies because the
relative contribution of each component shifts with stress. The high-strain and
Glen-Nye laws are grain-size insensitive single-component laws that produce straight
lines on a log–log plot, but with different exponents (``n = 3.5`` vs ``n = 3``).
**Right panel** — viscosity of the low-strain law as a function of grain size at
``\sigma_e = 0.3`` MPa. Coarser-grained ice is stiffer because the GSS mechanisms
slow down with increasing ``d``. The Glen-Nye law is grain-size insensitive.
=#

T_fan_K     = 268.15                                                  # −5 °C in K
σ_range_fan = 10 .^ range(log10(1e3), log10(1e6), length = 300)      # 1 kPa – 1 MPa

fan_ls_1mm = FanLowStrainFlowLaw(1e-3, T_fan_K)
fan_ls_5mm = FanLowStrainFlowLaw(5e-3, T_fan_K)
fan_hs     = FanHighStrainFlowLaw()
glen_nf    = GlenNyeFlowLaw()

A_fan_ls = rate_factor(T_fan_K, fan_ls_1mm.rate_factor)   # = 1.0; baked into creep
A_fan_hs = rate_factor(T_fan_K, fan_hs.rate_factor)
A_glen   = rate_factor(T_fan_K, glen_nf.rate_factor)

η_ls_1mm = viscosity(A_fan_ls, creep(σ_range_fan, fan_ls_1mm.creep), fan_ls_1mm)
η_ls_5mm = viscosity(A_fan_ls, creep(σ_range_fan, fan_ls_5mm.creep), fan_ls_5mm)
η_hs     = viscosity(A_fan_hs, creep(σ_range_fan, fan_hs.creep),     fan_hs)
η_glen   = viscosity(A_glen,   creep(σ_range_fan, glen_nf.creep),    glen_nf)

d_range  = 10 .^ range(log10(3e-4), log10(3e-2), length = 300)    # 0.3 mm – 30 mm
σ_fixed  = 0.3e6   # Pa

η_vs_d = [begin
              ls = FanLowStrainFlowLaw(d, T_fan_K)
              viscosity(rate_factor(T_fan_K, ls.rate_factor), creep(σ_fixed, ls.creep), ls)
          end for d in d_range]
η_glen_ref = viscosity(A_glen, creep(σ_fixed, glen_nf.creep), glen_nf)

fig_fan = Figure(size = (900, 420))

ax_σ = Axis(fig_fan[1, 1],
    xlabel = L"Effective stress $\sigma_e$ (kPa)",
    ylabel = L"Viscosity $\eta$ (Pa s)",
    xscale = log10,
    yscale = log10,
    title  = L"Viscosity vs stress at $T = {-5}\,$°C",
)
lines!(ax_σ, σ_range_fan ./ 1e3, η_ls_1mm, label = L"Fan low-strain, $d = 1$ mm")
lines!(ax_σ, σ_range_fan ./ 1e3, η_ls_5mm, label = L"Fan low-strain, $d = 5$ mm",    linestyle = :dash)
lines!(ax_σ, σ_range_fan ./ 1e3, η_hs,     label = "Fan high-strain (1-component)",  linestyle = :dashdot)
lines!(ax_σ, σ_range_fan ./ 1e3, η_glen,   label = "Glen-Nye",                       linestyle = :dot, color = :gray)
axislegend(ax_σ, position = :lb, labelsize = 11)

ax_d = Axis(fig_fan[1, 2],
    xlabel = L"Grain size $d$ (mm)",
    ylabel = L"Viscosity $\eta$ (Pa s)",
    xscale = log10,
    yscale = log10,
    title  = L"Grain-size sensitivity at $T = {-5}\,$°C, $\sigma_e = 0.3$ MPa",
)
lines!(ax_d, d_range .* 1e3, η_vs_d,    label = "Fan low-strain (3-component)")
hlines!(ax_d, [η_glen_ref],              label = "Glen-Nye (grain-size insensitive)", linestyle = :dot, color = :gray)
axislegend(ax_d, position = :lt, labelsize = 11)

fig_fan

#=
In [`AbstractFlowLaw`](@ref), we show convenience constructors, other options, as well as how to implement your own flow law.

## [Anisotropy](@id anisotropy)

Ice being a polycrystalline material, its deformation can be influenced by the orientation of its crystal fabric and the presence of impurities. These factors can lead to anisotropic behavior, where the ice deforms more easily in certain directions than others. Pagos.jl implements various parameterizations of this effect via [`AbstractAnisotropy`](@ref) and allows users to implement their own models.

For instance, [`EnhancementFactorAnisotropy`](@ref) applies a constant enhancement factor ``E`` to the flow law, softening or hardening the ice uniformly in all directions. This value typically ranges from 0.1 to 10 in glaciological applications and can be used to account for the effects of crystal fabric or impurities on ice deformation without explicitly modeling the microstructural details. The modified flow law incorporating the enhancement factor can be expressed as:

```math
\begin{aligned}
A(T') \rightarrow E A(T')
\end{aligned}
```

In a continental ice-sheet model, different enhacement factors are typically applied to the shear, stream and shelf regions to capture the effects of fabric development and impurity content on ice deformation. For instance, a common choice is to use a higher enhancement factor (e.g., ``E = 2``) in stream and shelf regions where the ice is softer due to fabric development, and a lower enhancement factor (e.g., ``E = 1``) in regions that were not subject yet to a lot of deformation, as it is typically the case in interior, shear-dominated regions.
=#

shear_anisotropy = EnhancementFactorAnisotropy(1.0)
stream_anisotropy = EnhancementFactorAnisotropy(2.0)
shelf_anisotropy  = EnhancementFactorAnisotropy(3.0)

#=
However, this approach does not capture the directional dependence of anisotropy, as it applies the same enhancement factor regardless of the loading direction. More sophisticated models, such as [`CAFFEAnisotropy`](@ref), explicitly account for the fabric state and its evolution under deformation, allowing for a more realistic representation of anisotropic behavior.

[`SimpleCAFFEAnisotropy`](@ref) offers a scalar simplification of the CAFFE model: the enhancement factor ``E`` is expressed as a function of a single deformability scalar ``\mathcal{D} \in [0, 5/2]``, which encodes the local fabric state. The piecewise law reads:

```math
E(\mathcal{D}) = \begin{cases}
    E_\min + (1 - E_\min)\,\mathcal{D}^{t} & 0 \le \mathcal{D} \le 1 \\[4pt]
    \dfrac{4\mathcal{D}^2(E_\max - 1) + 25 - 4 E_\max}{21} & 1 < \mathcal{D} \le \tfrac{5}{2}
\end{cases}
```

where ``t = \tfrac{8}{21}(E_\max - 1)/(1 - E_\min)`` and the default parameters are ``E_\min = 0.1``, ``E_\max = 10``. ``\mathcal{D} = 0`` corresponds to a fully hardened (single-maximum) fabric, ``\mathcal{D} = 1`` to isotropic fabric, and ``\mathcal{D} = 5/2`` to a fully soft (girdle) fabric.
=#

simple_caffe = SimpleCAFFEAnisotropy()
D_range = range(0, stop = 5/2, length = 500)
E_simple = map(D -> enhancement_factor(D, simple_caffe), D_range)

fig_simple_caffe = Figure()
ax_sc = Axis(fig_simple_caffe[1, 1],
    xlabel = L"Deformability $\mathcal{D}$",
    ylabel = L"Enhancement factor $E$",
    title  = "SimpleCAFFEAnisotropy enhancement factor",
)
lines!(ax_sc, collect(D_range), E_simple)
vlines!(ax_sc, [1.0], linestyle = :dash, color = :gray, label = "Isotropic fabric")
axislegend(ax_sc, position = :lt)
fig_simple_caffe

#=
## [Pressure melting point](@id melting_point)

For some computations, it is necessary to determine the temperature relative to the pressure melting point, ``T'``. The pressure melting point decreases with increasing pressure, and can be approximated using a linear relation:
=#

lpmp = LinearPressureMeltingPoint()
p = range(0, stop = 50f6, length = 1000)
Tm = map(x -> pressure_melting_point(x, lpmp), p)
fig_tm = plot_melting_point(p, Tm)

#=
In [`AbstractPressureMeltingPoint`](@ref), we show other options, as well as how to implement your own pressure melting point models.
=#

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
[`GlenNyeCreep`](@ref) is the standard power law ``f = \sigma_e^{n-1}`` (slope 2 for
``n = 3``).
[`RegularizedGlenNyeCreep`](@ref) adds a floor ``\sigma_0^{n-1}`` that prevents a
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
cf_glen    = creep(σ_range, GlenNyeCreep())
cf_reg     = creep(σ_range, RegularizedGlenNyeCreep(σ_0 = 1e3))
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
In [`AbstractFlowLaw`](@ref), we show convenience constructors, other options, as well as how to implement your own flow law.

## [Enhancement factor](@id enhancement_factor)

In glaciology, an enhancement factor ``E`` is often introduced to account for deviations from the standard flow law due to factors such as impurities, crystal orientation, or other microstructural effects. The modified flow law incorporating the enhancement factor can be expressed as:

```math
\begin{aligned}
A(T') \rightarrow E A(T')
\end{aligned}
```

The enhancement factor ``E`` is typically a dimensionless quantity greater than 1, indicating that the ice deforms more easily than predicted by the standard flow law. It can vary depending on the specific conditions and characteristics of the ice being studied.

In Pagos.jl, the enhancement factor is set through the `E` field of [`GlenNyeCreep`](@ref) and [`RegularizedGlenNyeCreep`](@ref). Values ``E > 1`` soften the ice (e.g. anisotropic crystal fabric, elevated impurity content), while ``E < 1`` harden it:
=#

glennye_soft = GlenNyeCreep(E = 3.0)
glennye_hard = GlenNyeCreep(E = 0.5)

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

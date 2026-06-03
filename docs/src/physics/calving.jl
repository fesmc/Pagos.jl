#=

# [Calving](@id calving)

Calving laws describe the process by which icebergs detach from an ice-sheet margin.
All laws are subtypes of [`AbstractCalving`](@ref) and are evaluated pointwise via
[`calving_rate`](@ref). The returned value is the local calving rate
``\dot{c}\,[\mathrm{m\,yr^{-1}}]``, with the convention that ``\dot{c} < 0`` removes ice.

Laws fall into three families depending on their physical driver:

| Family | Laws |
|--------|------|
| Ice thickness | [`PrescribedCalving`](@ref), [`RelaxedCalving`](@ref), [`ThicknessCalving`](@ref), [`FlotationCalving`](@ref) |
| Marine cliff geometry | [`PollardDeContoCalving`](@ref), [`CrawfordCalving`](@ref) |
| Flow divergence | [`EigenCalving`](@ref), [`LevermannCalving`](@ref), [`LipscombCalving`](@ref) |

The simplest default is [`PrescribedCalving`](@ref): the default `rate = 0` effectively
disables calving, providing a safe starting point when the calving boundary is not the
primary interest:

=#

using Pagos, CairoMakie

no_calving    = PrescribedCalving()               # ċ = 0 everywhere — calving disabled
const_calving = PrescribedCalving(rate = 10.0)    # uniform 10 m yr⁻¹

#=

For simulations where the calving front matters, [`ThicknessCalving`](@ref) is a
widely-used first choice. It removes ice in excess of a critical thickness
``H_\mathrm{crit}`` on a relaxation timescale ``\tau``, and is simple to tune from
observations or sensitivity studies:

=#

thickness_calving = ThicknessCalving()   # H_crit = 100 m, τ = 1 yr (defaults)

#=

Other laws depend on marine cliff geometry or flow divergence and require additional
inputs (surface elevation, bed elevation, or strain rates). See the comparison figure
below and the [`AbstractCalving`](@ref) API reference for all available subtypes.

!!! note "Laws requiring additional context"
    [`LipscombCalving`](@ref) and [`LevermannCalving`](@ref) access the local cell
    dimensions ``\Delta x, \Delta y`` as global variables in the scalar
    `calving_rate` method — use the in-place `calving_rate!` form with a full state.
    [`BassisCalving`](@ref) requires a constants struct carrying ``\rho_\mathrm{ice}``,
    ``\rho_\mathrm{seawater}``, and ``g``; use it via `calving_rate!` with a populated
    [`Params`](@ref).

## [Make your own calving law](@id custom_calving)

Any subtype of [`AbstractCalving`](@ref) can be used by implementing the scalar
[`calving_rate`](@ref) method. Below is an example of a law that drives calving
proportionally to ocean temperature above a threshold:

=#

@kwdef struct OceanTemperatureCalving{T} <: AbstractCalving
    T_threshold::T = 0.0      # °C — no calving below this temperature
    rate_per_degree::T = 5.0  # m yr⁻¹ °C⁻¹
    max_rate::T = Inf          # m yr⁻¹
end

function Pagos.calving_rate(T_ocean, calving::OceanTemperatureCalving)
    (; T_threshold, rate_per_degree, max_rate) = calving
    excess = max(T_ocean - T_threshold, 0.0)
    return -min(rate_per_degree * excess, max_rate)
end

ocean_calving = OceanTemperatureCalving(T_threshold = 1.0, rate_per_degree = 8.0)

#=

## Comparing all calving laws

The figure below visualises the three families. All panels show the calving magnitude
``|\dot{c}|\,[\mathrm{m\,yr^{-1}}]``; recall that the actual returned value is negative
(ice removed).

**Left — thickness-based laws** (sweep over ice thickness ``H``).
[`PrescribedCalving`](@ref) is a flat baseline. [`RelaxedCalving`](@ref) grows linearly
with ``H`` once the threshold is crossed, while [`ThicknessCalving`](@ref) grows with the
*excess* ``H - H_\mathrm{crit}`` — a softer response. [`FlotationCalving`](@ref) is
active only below flotation thickness (here ``H_\mathrm{float} = 400\,\mathrm{m}``) and
removes more ice as the calving front thickens toward the grounding threshold.

**Centre — marine cliff geometry laws** (sweep over subaerial cliff height
``H_s = z_\mathrm{srf} - z_\mathrm{sl}``). [`PollardDeContoCalving`](@ref) grows
linearly beyond its threshold. [`CrawfordCalving`](@ref) starts almost dormant but
overtakes Pollard–DeConto around ``H_s \approx 295\,\mathrm{m}`` and grows explosively —
the ``\alpha = 7.3`` exponent encodes the structural fragility of very tall ice cliffs.

**Right — strain-rate-based law** (sweep over ``\dot{\varepsilon}_1`` with fixed
``\dot{\varepsilon}_2``). [`EigenCalving`](@ref) produces a linear ramp in
``\dot{\varepsilon}_1`` whose slope is set by ``K \dot{\varepsilon}_2``: larger
transverse extension steepens the response proportionally.

=#

H_range = range(1.0, 600.0, length = 600)
H_float = 400.0

c_const = fill(abs(calving_rate(PrescribedCalving(rate = 10.0))), length(H_range))
c_relax = abs.([calving_rate(h, RelaxedCalving(H_critical = 200.0, timescale = 1.0)) for h in H_range])
c_thick = abs.([calving_rate(h, h, ThicknessCalving(H_critical = 200.0, timescale = 1.0)) for h in H_range])
c_float = abs.([calving_rate(h, h - H_float, FlotationCalving(timescale = 1.0)) for h in H_range])

H_s_range = range(80.0, 350.0, length = 500)
z_sl  = 0.0
z_bed = -1000.0

c_pollard  = abs.([calving_rate(h, z_sl, z_bed, PollardDeContoCalving(H_c = 100.0, timescale = 1.0)) for h in H_s_range])
c_crawford = abs.([calving_rate(h, z_sl, z_bed, CrawfordCalving()) for h in H_s_range])

ε₁_range = range(0.0, 1e-2, length = 500)

c_eig_lo = abs.([calving_rate(e, 1e-3, EigenCalving(K = 1e7)) for e in ε₁_range])
c_eig_mi = abs.([calving_rate(e, 5e-3, EigenCalving(K = 1e7)) for e in ε₁_range])
c_eig_hi = abs.([calving_rate(e, 1e-2, EigenCalving(K = 1e7)) for e in ε₁_range])

set_theme!(theme_latexfonts())
fig_calv = Figure(size = (1000, 360))

ax1 = Axis(fig_calv[1, 1],
    xlabel = L"Ice thickness $H$ (m)",
    ylabel = L"$|\dot{c}|$ (m yr$^{-1}$)",
    title  = "Thickness-based",
)
lines!(ax1, H_range, c_const, label = L"Prescribed (rate $= 10$ m yr$^{-1}$)")
lines!(ax1, H_range, c_relax, label = L"Relaxed ($H_\mathrm{crit} = 200$ m, $\tau = 1$ yr)")
lines!(ax1, H_range, c_thick, label = L"Thickness ($H_\mathrm{crit} = 200$ m, $\tau = 1$ yr)")
lines!(ax1, H_range, c_float, label = L"Flotation ($\tau = 1$ yr, $H_\mathrm{float} = 400$ m)", linestyle = :dash)
axislegend(ax1, position = :lt, labelsize = 9)

ax2 = Axis(fig_calv[1, 2],
    xlabel = L"Cliff height above sea level $H_s$ (m)",
    ylabel = L"$|\dot{c}|$ (m yr$^{-1}$)",
    title  = "Marine cliff geometry",
)
lines!(ax2, H_s_range, c_pollard,  label = L"Pollard-DeConto ($H_c = 100$ m, $\tau = 1$ yr)")
lines!(ax2, H_s_range, c_crawford, label = L"Crawford ($I = 1.9 \times 10^{-16}$, $\alpha = 7.3$)")
vlines!(ax2, [100.0, 135.0], color = :gray, linestyle = :dot)
ylims!(ax2, 0, 500)
axislegend(ax2, position = :lt, labelsize = 9)

ax3 = Axis(fig_calv[1, 3],
    xlabel = L"Extensional strain rate $\dot{\varepsilon}_1$ (yr$^{-1}$)",
    ylabel = L"$|\dot{c}|$ (m yr$^{-1}$)",
    title  = L"EigenCalving ($K = 10^7$ m yr)",
)
lines!(ax3, ε₁_range, c_eig_lo, label = L"$\dot{\varepsilon}_2 = 10^{-3}$ yr$^{-1}$")
lines!(ax3, ε₁_range, c_eig_mi, label = L"$\dot{\varepsilon}_2 = 5 \times 10^{-3}$ yr$^{-1}$")
lines!(ax3, ε₁_range, c_eig_hi, label = L"$\dot{\varepsilon}_2 = 10^{-2}$ yr$^{-1}$")
axislegend(ax3, position = :lt, labelsize = 9)

fig_calv

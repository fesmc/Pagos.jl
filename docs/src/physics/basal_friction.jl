#=

# [Basal friction](@id basal_friction)

Basal friction laws relate the basal shear stress ``\boldsymbol{\tau}_{\mathrm{b}}``
to the basal sliding velocity ``\mathbf{v}_{\mathrm{b}}``. All laws share the common
structure:

```math
\boldsymbol{\tau}_{\mathrm{b}} = -\beta(\mathbf{v}_{\mathrm{b}}) \, \mathbf{v}_{\mathrm{b}}
```

where ``\beta \geq 0`` is the basal friction coefficient computed by any subtype of
[`AbstractBasalBeta`](@ref). The laws differ in how ``\beta`` depends on the sliding
speed and on the bed strength ``c_{\mathrm{b}}``.

The recommended default is [`CoulombBasalBeta`](@ref), which saturates to the bed
strength ``c_{\mathrm{b}}`` at high sliding velocities — recovering the classical
Coulomb friction limit — while remaining smooth and well-behaved at low velocities:

=#

using Pagos, CairoMakie

coulomb_beta = CoulombBasalBeta()   # q = 0.2, v₀ = 100 m yr⁻¹ (defaults)

#=

The simplest alternative is [`ConstantBasalBeta`](@ref), which fixes ``\beta`` to a
uniform value and gives linear (viscous) sliding with no velocity dependence or
saturation. [`PseudoPlasticPowerBasalBeta`](@ref) offers a power-law regime without
the Coulomb plateau: it recovers perfectly plastic sliding at ``q = 0`` and linear
sliding at ``q = 1``.

## [Make your own basal friction law](@id custom_beta)

Any subtype of [`AbstractBasalBeta`](@ref) can be used by implementing the scalar
[`basal_beta`](@ref) method. The generic `basal_beta!` fallback then handles field
operations automatically. Below is an example of a law that is frictionless below a
threshold sliding speed and applies a constant ``\beta`` above it:

=#

@kwdef struct ThresholdBasalBeta{T} <: AbstractBasalBeta
    v_threshold::T = 10.0   # m yr⁻¹ — no friction below this speed
    β::T = 1e3              # Pa yr m⁻¹
end

function Pagos.basal_beta(_c_bed, v_basal, bb::ThresholdBasalBeta)
    return norm(v_basal) < bb.v_threshold ? zero(bb.β) : bb.β
end

threshold_beta = ThresholdBasalBeta()

#=

## Comparing all parameterizations

The figure below shows the normalized stress magnitude
``|\boldsymbol{\tau}_{\mathrm{b}}| / c_{\mathrm{b}}`` as a function of sliding speed for
the default reference velocity ``v_0 = 100\,\mathrm{m\,yr^{-1}}``.

**Left — power law** ([`PseudoPlasticPowerBasalBeta`](@ref)).
At ``q = 0`` the law is perfectly plastic: the stress equals ``c_{\mathrm{b}}``
everywhere, independent of velocity. Increasing ``q`` toward 1 transitions the response
toward linear sliding, where stress grows proportionally with speed with no upper bound.

**Centre — Coulomb, varying ``q``** ([`CoulombBasalBeta`](@ref)).
All curves saturate to ``c_{\mathrm{b}}`` at high velocities — the defining property of
Coulomb friction. The default ``q = 0.2`` gives a gentle, near-flat response even at
moderate speeds, while ``q = 1`` is the classical Coulomb form.

**Right — effect of ``v_{\mathrm{reg}}``** ([`CoulombBasalBeta`](@ref), ``q = 0.2``).
Without regularization (``v_{\mathrm{reg}} = 0``) the stress is strictly zero at rest.
Adding ``v_{\mathrm{reg}} > 0`` following [zoet_slip_2020](@citet) shifts the onset of
friction to finite velocity and smooths the near-zero behaviour, aiding numerical
convergence. The Coulomb plateau is unaffected.

The dotted vertical line marks ``v_0``, and the dashed horizontal line marks the
saturation stress ``c_{\mathrm{b}}``.

=#

v_range = range(0.0, 300.0, length = 500)   # m yr⁻¹
v_0  = 100.0                                 # default transition velocity
q_df = 0.2                                   # default Coulomb exponent

τ_power_q00 = ones(length(v_range))
τ_power_q05 = (v_range ./ v_0) .^ 0.5
τ_power_q08 = (v_range ./ v_0) .^ 0.8
τ_power_q10 = (v_range ./ v_0) .^ 1.0

τ_coulomb_q02 = (v_range ./ (v_range .+ v_0)) .^ 0.2
τ_coulomb_q05 = (v_range ./ (v_range .+ v_0)) .^ 0.5
τ_coulomb_q10 = (v_range ./ (v_range .+ v_0)) .^ 1.0

coulomb_norm(v, vr) = ((v + vr) / (v + vr + v_0)) ^ q_df
τ_coulomb_vreg0  = coulomb_norm.(v_range, 0.0)
τ_coulomb_vreg1  = coulomb_norm.(v_range, 1.0)
τ_coulomb_vreg5  = coulomb_norm.(v_range, 5.0)
τ_coulomb_vreg20 = coulomb_norm.(v_range, 20.0)

set_theme!(theme_latexfonts())
fig_bf = Figure(size = (1100, 380))

ax1 = Axis(fig_bf[1, 1],
    xlabel = L"Sliding speed $|v_b|$ (m yr$^{-1}$)",
    ylabel = L"$|\tau_b| / c_b$",
    title  = "PseudoPlasticPowerBasalBeta",
    limits = (0, 300, 0, 3),
)
lines!(ax1, v_range, τ_power_q00, label = L"$q = 0$ (plastic)",   linestyle = :dash)
lines!(ax1, v_range, τ_power_q05, label = L"$q = 0.5$")
lines!(ax1, v_range, τ_power_q08, label = L"$q = 0.8$ (default)")
lines!(ax1, v_range, τ_power_q10, label = L"$q = 1$ (linear)",    linestyle = :dot)
vlines!(ax1, [v_0], color = :gray, linestyle = :dot, label = L"$v_0 = 100$ m yr$^{-1}$")
axislegend(ax1, position = :lt, labelsize = 10)

ax2 = Axis(fig_bf[1, 2],
    xlabel = L"Sliding speed $|v_b|$ (m yr$^{-1}$)",
    ylabel = L"$|\tau_b| / c_b$",
    title  = L"CoulombBasalBeta ($v_\mathrm{reg} = 0$)",
    limits = (0, 200, 0, 1.2),
)
lines!(ax2, v_range, τ_coulomb_q02, label = L"$q = 0.2$ (default)")
lines!(ax2, v_range, τ_coulomb_q05, label = L"$q = 0.5$")
lines!(ax2, v_range, τ_coulomb_q10, label = L"$q = 1$ (standard Coulomb)", linestyle = :dot)
hlines!(ax2, [1.0], color = :gray, linestyle = :dash, label = L"Coulomb limit $c_b$")
vlines!(ax2, [v_0], color = :gray, linestyle = :dot, label = L"$v_0 = 100$ m yr$^{-1}$")
axislegend(ax2, position = :rb, labelsize = 10)

ax3 = Axis(fig_bf[1, 3],
    xlabel = L"Sliding speed $|v_b|$ (m yr$^{-1}$)",
    ylabel = L"$|\tau_b| / c_b$",
    title  = "CoulombBasalBeta",
    limits = (0, 200, 0, 1.2),
)
lines!(ax3, v_range, τ_coulomb_vreg0,  label = L"$v_\mathrm{reg} = 0$ (default)")
lines!(ax3, v_range, τ_coulomb_vreg1,  label = L"$v_\mathrm{reg} = 1$ m yr$^{-1}$")
lines!(ax3, v_range, τ_coulomb_vreg5,  label = L"$v_\mathrm{reg} = 5$ m yr$^{-1}$")
lines!(ax3, v_range, τ_coulomb_vreg20, label = L"$v_\mathrm{reg} = 20$ m yr$^{-1}$")
hlines!(ax3, [1.0], color = :gray, linestyle = :dash, label = L"Coulomb limit $c_b$")
vlines!(ax3, [v_0], color = :gray, linestyle = :dot, label = L"$v_0 = 100$ m yr$^{-1}$")
axislegend(ax3, position = :rb, labelsize = 10)

fig_bf

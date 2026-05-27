#=

# [Topography](@id topography)

## [Sigma transform](@id sigma_transform)

Pagos.jl uses terrain-following (sigma) coordinates for the vertical dimension.
The sigma coordinate ``\sigma \in [0, 1]`` maps each ice column onto the unit interval
with ``\sigma = 0`` at the bed and ``\sigma = 1`` at the surface. A given sigma level
corresponds to the physical elevation:

```math
z = z_{\mathrm{b}} + \sigma \, H
```

where ``z_{\mathrm{b}}`` is the bed elevation and ``H`` the local ice thickness.
Any subtype of [`AbstractSigmaTransform`](@ref) specifies how the ``n`` layer midpoints
(``\zeta_{\mathrm{aa}}``) and ``n+1`` layer interfaces (``\zeta_{\mathrm{ac}}``) are
distributed within ``[0, 1]``, allowing resolution to be concentrated near the bed
or the surface.

The recommended default is [`CorrectedVerticalLayering`](@ref) with
[`PowerSigmaTransform`](@ref) at ``p = 1`` (uniform spacing). Unlike
[`VerticalLayering`](@ref), it places interfaces first and derives midpoints as their
arithmetic means, so each cell centre is the exact geometric midpoint of its bounding
faces — the natural choice for finite-volume schemes:

=#

using Pagos, CairoMakie

T = Float32
default_layering = CorrectedVerticalLayering(T, PowerSigmaTransform(10, 1))

#=

Increasing the exponent ``p > 1`` concentrates layers near the bed, which is useful
when higher resolution is needed to resolve basal processes or temperature gradients:

=#

quadratic_layering = CorrectedVerticalLayering(T, PowerSigmaTransform(10, 2))

#=

## [Make your own sigma transform](@id custom_sigma)

Any subtype of [`AbstractSigmaTransform`](@ref) can be used with [`VerticalLayering`](@ref)
by extending `Pagos.get_ζ_aa`. Below is an exponential transform that concentrates layers
more aggressively near the bed than a power law:

=#

struct ExponentialSigmaTransform <: AbstractSigmaTransform
    n::Int
    scale::Float64   # larger → stronger near-bed concentration
end
ExponentialSigmaTransform(; n = 10, scale = 3.0) = ExponentialSigmaTransform(n, scale)

function Pagos.get_ζ_aa(T, transform::ExponentialSigmaTransform)
    (; n, scale) = transform
    t = range(0.0, 1.0, length = n)
    ζ = (exp.(scale .* t) .- 1) ./ (exp(scale) - 1)
    return T.(ζ)
end

exp_layering = VerticalLayering(T, ExponentialSigmaTransform(n = 10, scale = 3.0))

#=

## Comparing vertical layerings

[`VerticalLayering`](@ref) places midpoints first (via the transform) and derives
interfaces as arithmetic means of adjacent midpoints. [`CorrectedVerticalLayering`](@ref)
inverts this: interfaces are placed first by the power law and midpoints are set as their
means. The difference is negligible for ``p = 1`` and becomes visible for ``p = 2``, where
the corrected variant produces a more regular distribution.

The figure below shows ``\zeta_{\mathrm{aa}}`` and ``\zeta_{\mathrm{ac}}`` for both
variants and two exponents.

=#

labels = [
    L"$\zeta_\mathrm{aa}, \: p = 1$",
    L"$\zeta_\mathrm{ac}, \: p = 1$",
    L"$\zeta_\mathrm{aa}, \: p = 2$",
    L"$\zeta_\mathrm{ac}, \: p = 2$",
]

linear_v    = VerticalLayering(T, PowerSigmaTransform(10, 1))
quadratic_v = VerticalLayering(T, PowerSigmaTransform(10, 2))
linear_vc    = CorrectedVerticalLayering(T, PowerSigmaTransform(10, 1))
quadratic_vc = CorrectedVerticalLayering(T, PowerSigmaTransform(10, 2))

set_theme!(theme_latexfonts())
fig = Figure(size = (800, 380))
ax1 = Axis(fig[1, 1], title = "VerticalLayering",          xlabel = "Layer index", ylabel = L"$\sigma$")
ax2 = Axis(fig[1, 2], title = "CorrectedVerticalLayering", xlabel = "Layer index", ylabel = L"$\sigma$")

for (ax, lin, quad) in [(ax1, linear_v, quadratic_v), (ax2, linear_vc, quadratic_vc)]
    scatterlines!(ax, lin.ζ_aa,  label = labels[1])
    scatterlines!(ax, lin.ζ_ac,  label = labels[2])
    scatterlines!(ax, quad.ζ_aa, label = labels[3], linestyle = :dash)
    scatterlines!(ax, quad.ζ_ac, label = labels[4], linestyle = :dash)
    axislegend(ax, position = :rb, labelsize = 10)
end

fig

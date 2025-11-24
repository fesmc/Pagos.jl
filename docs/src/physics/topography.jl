#=

# [Topography](@ref topography_examples)

## Sigma transform

The sigma transform is a common technique used in glaciological modeling to map the physical vertical coordinate ``z`` to a normalized vertical coordinate ``\sigma`` that ranges from 0 at the ice surface to 1 at the ice base. This transformation simplifies the representation of the ice geometry and allows for more efficient numerical computations.

=#

using Pagos, CairoMakie

T = Float32
linear_layering = VerticalLayering(T, PowerSigmaTransform(10, 1))
quadratic_layering = VerticalLayering(T, PowerSigmaTransform(10, 2))
labels = [
    L"$\zeta_\mathrm{aa}, \: p = 1$",
    L"$\zeta_\mathrm{ac}, \: p = 1$",
    L"$\zeta_\mathrm{aa}, \: p = 2$",
    L"$\zeta_\mathrm{ac}, \: p = 2$",
]
fig, ax, _ = scatterlines(linear_layering.ζ_aa, label = labels[1])
scatterlines!(ax, linear_layering.ζ_ac, label = labels[2])
scatterlines!(ax, quadratic_layering.ζ_aa, label = labels[3])
scatterlines!(ax, quadratic_layering.ζ_ac, label = labels[4])
axislegend(ax, position = :rb)
fig

linear_layering = CorrectedVerticalLayering(T, PowerSigmaTransform(10, 1))
quadratic_layering = CorrectedVerticalLayering(T, PowerSigmaTransform(10, 2))
fig, ax, _ = scatterlines(linear_layering.ζ_aa, label = labels[1])
scatterlines!(ax, linear_layering.ζ_ac, label = labels[2])
scatterlines!(ax, quadratic_layering.ζ_aa, label = labels[3])
scatterlines!(ax, quadratic_layering.ζ_ac, label = labels[4])
axislegend(ax, position = :rb)
fig
"""
    sigma(z, b, H)

Return the terrain-following sigma coordinate of a point at physical elevation `z`,
given bed elevation `b` and ice thickness `H`.

The sigma coordinate is defined as

```math
\\sigma = \\frac{z - b}{H}
```

so that ``\\sigma = 0`` at the bed and ``\\sigma = 1`` at the ice surface.
"""
sigma(z, b, H) = (z - b) / H

"""
    exponential_vertical_layers(n)

Return a vector of `n` sigma-coordinate layer boundaries with exponential spacing.

The layer boundaries are mapped from the uniform grid `[1/n, 2/n, ..., 1]` via

```math
\\sigma_l = \\frac{e^{l/n} - 1}{e - 1}, \\quad l = 1, \\ldots, n
```

so that the spacing increases toward the surface (``\\sigma = 1``), giving finer
resolution near the bed where shear deformation is largest.
"""
function exponential_vertical_layers(n)
    dz = 1 / n
    return (exp.(dz:dz:1) .- exp(0)) ./ (exp(1) - exp(0))
end
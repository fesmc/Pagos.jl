"""
$(TYPEDSIGNATURES)

Abstract type for terrain-following vertical coordinate transforms.

The sigma (``\\sigma``) coordinate maps the ice column onto the unit interval
``[0, 1]``, with ``\\sigma = 0`` at the bed and ``\\sigma = 1`` at the surface.
For a column of thickness ``H`` above a bed at elevation ``z_{\\mathrm{b}}``,
the physical elevation ``z`` corresponding to a given ``\\sigma`` level is:

```math
z = z_{\\mathrm{b}} + \\sigma \\, H
```

Subtypes specify how layer interfaces (``\\zeta_{\\mathrm{ac}}``) and layer
midpoints (``\\zeta_{\\mathrm{aa}}``) are distributed within ``[0, 1]``,
allowing resolution to be concentrated near the bed or the surface.

# Available subtypes
 - [`PowerSigmaTransform`](@ref)
 - [`ArctanSigmaTransform`](@ref)
"""
abstract type AbstractSigmaTransform end

"""
$(TYPEDSIGNATURES)

Sigma-level distribution based on a power law.

Layer midpoints follow:

```math
\\zeta_i = \\left(\\frac{i - 1}{n - 1}\\right)^{p}, \\quad i = 1, \\ldots, n
```

where ``p`` is the `exponent`. Values ``p > 1`` concentrate layers near the bed
(small ``\\zeta``); ``p = 1`` gives uniform spacing (see `LinearSigmaTransform`);
``p < 1`` concentrates layers near the surface.

# Fields
 - `n::Int`: number of vertical layers.
 - `exponent::T`: power-law exponent ``p``.

# Convenience constructors
 - `LinearSigmaTransform(T, n)` — uniform spacing (``p = 1``).
 - `QuadraticSigmaTransform(T, n)` — quadratic clustering near the bed (``p = 2``).
"""
struct PowerSigmaTransform{T} <: AbstractSigmaTransform
    n::Int
    exponent::T
end

"""
$(TYPEDSIGNATURES)

A convience constructor for `PowerSigmaTransform` with a linear distribution of layers (``p = 1``).
"""
LinearSigmaTransform(T, n) = PowerSigmaTransform{T}(n, 1)

"""
$(TYPEDSIGNATURES)

A convience constructor for `PowerSigmaTransform` with a quadratic distribution of layers (``p = 2``).
"""
QuadraticSigmaTransform(T, n) = PowerSigmaTransform{T}(n, 2)

"""
$(TYPEDSIGNATURES)

Sigma-level distribution based on an arctangent stretching.

Layers are clustered near the bed by mapping the unit interval through an
arctangent function scaled by `stretch_factor`. Larger values of
`stretch_factor` increase the concentration of layers near ``\\sigma = 0``.

# Fields
 - `n::Int`: number of vertical layers.
 - `stretch_factor::T`: controls the degree of near-bed clustering.
"""
struct ArctanSigmaTransform{T} <: AbstractSigmaTransform
    n::Int
    stretch_factor::T
end

function get_ζ_aa(T, transform::PowerSigmaTransform)
    (; n, exponent) = transform
    ζ_aa = range(0.0, stop = 1.0, length = n) .^ exponent
    return T.(ζ_aa)
end

function get_ζ_ac(ζ_aa)
    n = length(ζ_aa)
    ζ_ac = zeros(eltype(ζ_aa), n + 1)
    for i in 2:n
        ζ_ac[i] = 0.5 * (ζ_aa[i - 1] + ζ_aa[i])
    end
    ζ_ac[n+1] = 1
    return ζ_ac
end

"""
$(TYPEDSIGNATURES)

Vertical layer structure derived from an [`AbstractSigmaTransform`](@ref).

Layer midpoints (``\\zeta_{\\mathrm{aa}}``) are placed first using the transform,
and layer interfaces (``\\zeta_{\\mathrm{ac}}``) are then computed as arithmetic
means of adjacent midpoints, with ``\\zeta_{\\mathrm{ac},0} = 0`` and
``\\zeta_{\\mathrm{ac},n} = 1`` fixed.

# Fields
 - `transform::S`: the underlying [`AbstractSigmaTransform`](@ref).
 - `ζ_aa::Vector{T}`: sigma coordinates at layer midpoints (cell centres), length ``n``.
 - `ζ_ac::Vector{T}`: sigma coordinates at layer interfaces (cell faces), length ``n + 1``.

See also [`CorrectedVerticalLayering`](@ref) for the face-first variant.
"""
struct VerticalLayering{T, S}
    transform::S
    ζ_aa::Vector{T}
    ζ_ac::Vector{T}
end

function VerticalLayering(T, transform)
    ζ_aa = get_ζ_aa(T, transform)
    ζ_ac = get_ζ_ac(ζ_aa)
    return VerticalLayering(transform, ζ_aa, ζ_ac)
end

"""
$(TYPEDSIGNATURES)

Vertical layer structure in which cell faces are placed first by the transform
and cell centres are derived as their arithmetic means.

Unlike [`VerticalLayering`](@ref), which distributes midpoints and derives faces,
this variant applies the power law directly to the ``n + 1`` face positions:

```math
\\zeta_{\\mathrm{ac},i} = \\left(\\frac{i}{n}\\right)^{p}, \\quad i = 0, \\ldots, n
```

and then sets each midpoint as:

```math
\\zeta_{\\mathrm{aa},i} = \\frac{\\zeta_{\\mathrm{ac},i} + \\zeta_{\\mathrm{ac},i+1}}{2}
```

This guarantees that each cell centre lies at the true geometric midpoint of
its bounding faces in sigma space, which is the natural choice for
finite-volume discretizations.

# Fields
 - `n::Int`: number of vertical layers.
 - `ζ_aa::Vector{T}`: sigma coordinates at layer midpoints (cell centres), length ``n``.
 - `ζ_ac::Vector{T}`: sigma coordinates at layer interfaces (cell faces), length ``n + 1``.
"""
struct CorrectedVerticalLayering{T}
    n::Int
    ζ_aa::Vector{T}
    ζ_ac::Vector{T}
end

function CorrectedVerticalLayering(T, transform)
    ζ_ac = range(0.0, stop = 1.0, length = transform.n + 1) .^ transform.exponent
    ζ_aa = zeros(T, transform.n)
    for i in 1:transform.n
        ζ_aa[i] = 0.5 * (ζ_ac[i] + ζ_ac[i + 1])
    end
    return CorrectedVerticalLayering(transform.n, ζ_aa, T.(ζ_ac))
end

"""
$(TYPEDSIGNATURES)

Map a physical elevation `z` to the sigma coordinate ``\\sigma \\in [0, 1]``:

```math
\\sigma = \\frac{z - z_{\\mathrm{b}}}{H}
```

where `b` is the bed elevation ``z_{\\mathrm{b}}`` and `H` is the local ice thickness.
Returns 0 at the bed and 1 at the ice surface.
"""
sigma(z, b, H) = (z - b) / H

# function sigma_transform(transform::LinearSigmaTransform)
#     return range(0.0, 1.0; length=transform.n)
# end


# function sigma_transform(transform::ExponentialSigmaTransform)
#     dz = 1 / transform.n
#     return (exp.(dz:dz:1) .- exp(0)) ./ (exp(1) - exp(0))
# end
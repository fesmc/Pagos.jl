"""
    NormalizedGaussLegendre{T<:AbstractFloat}

A type that stores the nodes and weights of a Gauss-Legendre quadrature rule
on the interval [-1, 1].

# Fields
- `n::Int`: Number of nodes.
- `x::Vector{T}`: Nodes.

# Examples
```
ngl = NormalizedGaussLegendre(3)
NormalizedGaussLegendre{Float64}(3, [-0.7745966692414834, 0.0, 0.7745966692414834], [0.5555555555555556, 0.8888888888888888, 0.5555555555555556])
```
"""
struct NormalizedGaussLegendre{T<:AbstractFloat}
    n::Int
    x::Vector{T}
    w::Vector{T}
end

function NormalizedGaussLegendre(n::Int)
    x, w = gausslegendre(n)
    return NormalizedGaussLegendre(n, x, w)
end

"""
    scaled_gausslegendre(ngl, a, b)

Return the nodes and weights of a Gauss-Legendre quadrature rule on the interval
[a, b] given the nodes and weights of a Gauss-Legendre quadrature rule on the
interval [-1, 1].

# Arguments
- `ngl::NormalizedGaussLegendre{T}`: Normalized Gauss-Legendre quadrature rule.
- `a::T`: Lower bound of the interval.
- `b::T`: Upper bound of the interval.

# Examples
```
ngl = NormalizedGaussLegendre(3)
x, w = scaled_gausslegendre(ngl, 0.0, 1.0)
f(x) = x^2
sum(f.(x) .* w) ≈ 1/3
```
"""
function scaled_gausslegendre(ngl::NormalizedGaussLegendre{T}, a::T, b::T) where {T<:AbstractFloat}
    x = ngl.x .* (b - a) / 2 .+ (a + b) / 2
    w = ngl.w * (b - a) / 2
    return x, w
end
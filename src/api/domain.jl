"""
    Domain(T, lx, ly, dx, dy)

Define the domain of the ice sheet model, which contains:
- `nx::Int`: the number of grid points in the x-direction.
- `ny::Int`: the number of grid points in the y-direction.
- `dx::T`: the grid spacing in the x-direction.
- `dy::T`: the grid spacing in the y-direction.
- `x::Vector{T}`: the grid points in the x-direction.
- `y::Vector{T}`: the grid points in the y-direction.
- `lx::T`: the length of the domain in the x-direction.
- `ly::T`: the length of the domain in the y-direction.
- `null::Matrix{T}`: a matrix of zeros of size `nx` by `ny`.
- `X::Matrix{T}`: a matrix of x-coordinates of size `nx` by `ny`.
- `Y::Matrix{T}`: a matrix of y-coordinates of size `nx` by `ny`.

Example:

```julia
domain = Domain(Float64, 6000, 6000, 16, 16)
```
"""
struct Domain{T<:AbstractFloat}
    nx::Int
    ny::Int
    dx::T
    dy::T
    x::Vector{T}
    y::Vector{T}
    lx::T
    ly::T
    null::Matrix{T}
    X::Matrix{T}
    Y::Matrix{T}
end

function Domain(T, lx, ly, dx, dy)
    return Domain(T(lx), T(ly), T(dx), T(dy))
end

function Domain(lx::T, ly::T, dx::T, dy::T) where {T<:AbstractFloat}
    x = collect(range(T(0), step=dx, stop=lx)) .- lx/2
    y = collect(range(T(0), step=dy, stop=ly)) .- ly/2
    nx = length(x)
    ny = length(y)
    null = zeros(T, nx, ny)
    X = x * ones(T, ny)'
    Y = ones(T, nx) * y'
    return Domain(nx, ny, dx, dy, x, y, lx, ly, null, X, Y)
end

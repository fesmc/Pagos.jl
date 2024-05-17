"""
    Options{T<:AbstractFloat}

Determine the options for the ice sheet model, which includes:
 - `theta_v::T = T(0.6)`: the velocity relaxation parameter of the PT method.
 - `theta_mu::T = T(0.1)`: the viscosity relaxation parameter of the PT method.
 - `reltol::T = T(1e-3)`: the relative tolerance of the PT method.
 - `maxiter::Int = 100`: the maximum number of iterations of the PT method.

To initialize with default values:
```julia
options = Options{Float64}()
```

To initialize with custom values:
```julia
options = Options{Float64}(theta_v=0.5, theta_mu=0.1, reltol=1e-4, maxiter=200)
```
"""
@kwdef struct Options{T<:AbstractFloat}
    theta_v::T = T(0.6)
    theta_mu::T = T(0.1)
    reltol::T = T(1e-3)
    maxiter::Int = 100
    debug::Bool = false
end
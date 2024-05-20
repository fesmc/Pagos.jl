"""
    Params{T<:AbstractFloat}

Struct containing the parameters of the ice sheet model, which contains:
- `rho_ice::T`: the density of ice.
- `muB::T`: the bulk viscosity of ice.
- `g::T`: the gravitational acceleration.

Example using default values:
```julia
params = Params{Float64}()
```

Example using user-defined values:
```julia
params = Params{Float64}(910.0, 0.5, 9.81)
```
"""
@kwdef struct Params{T<:AbstractFloat}
    rho_ice::T = T(910.0)       # ice density
    muB::T = T(1e2)             # bulk viscosity of ice
    g::T = T(9.81)              # gravitational acceleration
end
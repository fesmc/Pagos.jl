"""
    Params{T<:AbstractFloat}

Struct containing the parameters of the ice sheet model, which contains:
- `ndim1::Int`: the dimension coefficient for 1D problems.
- `ndim2::Int`: the dimension coefficient for 2D problems.
- `ndim3::Int`: the dimension coefficient for 3D problems.
- `rho_ice::T`: the density of ice.
- `min_bulk_viscosity_ice::T`: the minimal bulk viscosity of ice.
- `muB::T`: the bulk viscosity of ice.
- `g::T`: the gravitational acceleration.

Example using default values:
```julia
params = Params{Float64}()
```

Example using user-defined values:
```julia
params = Params{Float64}(rho_ice = 910.0)
```
"""
@kwdef struct Params{T<:AbstractFloat}
    ndim1::T = 2.1
    ndim2::T = 4.1
    ndim3::T = 6.1
    rho_ice::T = 910.0               # ice density
    min_bulk_viscosity_ice::T = 0.5  # minimal bulk viscosity of ice
    muB::T = 1e2                     # (old) bulk viscosity of ice
    g::T = 9.81                      # gravitational acceleration
end

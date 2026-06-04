"""
    IceSheet{T}(state::State{T}, domain::Domain{T}, params::Params{T}, options::Options{T})

Struct containing the ice sheet model, which contains:
- `state::State{T}`: the [`State`](@ref) of the ice sheet model.
- `domain::Domain{T}`: the [`Domain`](@ref) of the ice sheet model.
- `params::Params{T}`: the [`Params`](@ref) of the ice sheet model.
- `options::Options{T}`: the [`Options`](@ref) of the ice sheet model.

Example:

```julia
T = Float64
domain = Domain(T, 6000, 6000, 16, 16)
state = State(domain)
params = Params{T}()
options = Options{T}()
icesheet = IceSheet(state, domain, params, options)
```
"""
struct IceSheet{T<:AbstractFloat}
    state::State{T}
    domain::Domain{T}
    params::Params{T}
    options::Options{T}
end
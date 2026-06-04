"""
    IceSheet(state, domain, params, options)

Struct containing the ice sheet model:
- `state::State{M}`: the [`State`](@ref) of the ice sheet model.
- `domain::Domain{T, V, M}`: the [`Domain`](@ref) of the ice sheet model.
- `params::Params{T}`: the [`Params`](@ref) of the ice sheet model.
- `options::Options{T}`: the [`Options`](@ref) of the ice sheet model.

CPU example:

```julia
T      = Float64
domain = Domain(T, 6000.0, 6000.0, 16.0, 16.0)
state  = State(domain)
params = Params{T}()
options = Options{T}()
icesheet = IceSheet(state, domain, params, options)
```

GPU example (requires a CUDA-capable device):

```julia
using CUDA, KernelAbstractions
domain_gpu = Domain(CUDABackend(), Float32, 6000.0, 6000.0, 16.0, 16.0)
state_gpu  = State(domain_gpu)
icesheet_gpu = IceSheet(state_gpu, domain_gpu, Params{Float32}(), Options{Float32}())
```
"""
struct IceSheet{S, D, P, O}
    state::S    # <:State
    domain::D   # <:Domain
    params::P   # <:Params
    options::O  # <:Options
end
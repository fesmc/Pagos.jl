"""
$(TYPEDSIGNATURES)

Define the domain of the ice sheet model, which contains:
- `nx`: the number of grid points in the x-direction.
- `ny`: the number of grid points in the y-direction.
- `dx`: the grid spacing in the x-direction.
- `dy`: the grid spacing in the y-direction.
- `x`: the grid points in the x-direction.
- `y`: the grid points in the y-direction.
- `lx`: the length of the domain in the x-direction.
- `ly`: the length of the domain in the y-direction.
- `null`: a matrix of zeros of size `nx` by `ny`.
- `X`: a matrix of x-coordinates of size `nx` by `ny`.
- `Y`: a matrix of y-coordinates of size `nx` by `ny`.

The two-argument form allocates CPU (`Array`) arrays.  Pass a
`KernelAbstractions.Backend` as the first argument to allocate on a different
device:

```julia
domain     = Domain(Float64, 6000.0, 6000.0, 16.0, 16.0)          # CPU
domain_gpu = Domain(CUDABackend(), Float64, 6000.0, 6000.0, 16.0, 16.0)  # GPU
```
"""
struct Domain{T, V, M}
    nx::Int
    ny::Int
    dx::T
    dy::T
    x::V
    y::V
    lx::T
    ly::T
    null::M
    X::M
    Y::M
end

function Domain(T::Type{<:AbstractFloat}, lx, ly, dx, dy)
    return Domain(CPU(), T, lx, ly, dx, dy)
end

function Domain(backend::Backend, T::Type{<:AbstractFloat}, lx, ly, dx, dy)
    x_cpu = collect(range(T(0), step = T(dx), stop = T(lx))) .- T(lx) / 2
    y_cpu = collect(range(T(0), step = T(dy), stop = T(ly))) .- T(ly) / 2
    nx    = length(x_cpu)
    ny    = length(y_cpu)
    x     = KernelAbstractions.adapt(backend, x_cpu)
    y     = KernelAbstractions.adapt(backend, y_cpu)
    null  = KernelAbstractions.zeros(backend, T, nx, ny)
    X     = KernelAbstractions.adapt(backend, x_cpu * ones(T, ny)')
    Y     = KernelAbstractions.adapt(backend, ones(T, nx) * y_cpu')
    return Domain(nx, ny, T(dx), T(dy), x, y, T(lx), T(ly), null, X, Y)
end

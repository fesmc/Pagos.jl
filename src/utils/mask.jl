"""
$(TYPEDSIGNATURES)

Pre-computed indices of cells where a boolean mask is `true`. Avoids branch divergence
in GPU kernels by replacing conditional iteration with direct indexing into a compact
list. Rebuild when the mask topology changes (e.g. after a calving event).

Move to a device with `adapt(backend, map)` from Adapt.jl; `get_backend(map.indices)`
then dispatches [`apply!`](@ref) to the correct kernel automatically.
"""
mutable struct ActiveCellsMap{I}
    indices::I
end

function ActiveCellsMap(mask::AbstractArray{Bool})
    return ActiveCellsMap(findall(mask))
end

function active_indices!(acm::ActiveCellsMap, mask::AbstractArray{Bool})
    acm.indices = findall(mask)
    return
end

@kernel function _apply_kernel!(z, f, active_map, arrays, scalars)
    n = @index(Global, Linear)
    I = active_map[n]
    @inbounds z[I] = f(ntuple(k -> arrays[k][I], Val(length(arrays)))..., scalars...)
end

"""
$(TYPEDSIGNATURES)

Apply `f` over all active indices in `map`, writing results into `z`:

```julia
for I in map
    z[I] = f(x[I], y[I], ..., p, q, ...)
end
```

Pass arrays to be indexed as a `Tuple`; any trailing arguments are forwarded to `f`
unchanged and can be used for multiple dispatch:

```julia
apply!(f, map, z, (x, y), (p,))   # z[I] = f(x[I], y[I], p)
```

The backend (CPU or GPU) is inferred from `map.indices`, so moving the index array to
the device is the only requirement for GPU dispatch.

`f` and all scalar/param arguments must be plain named functions or `isbits` values —
closures capture heap-allocated state and are not GPU-compatible.
"""
function apply!(f, map::ActiveCellsMap, z, arrays::Tuple, scalars::Tuple=(); workgroupsize=256)
    backend = get_backend(map.indices)
    kernel! = _apply_kernel!(backend, workgroupsize)
    kernel!(z, f, map.indices, arrays, scalars; ndrange=length(map.indices))
    KernelAbstractions.synchronize(backend)
    return
end

# Convenience method: when all inputs are arrays and there are no scalars.
function apply!(f, map::ActiveCellsMap, z, arrays::AbstractArray...; workgroupsize=256)
    apply!(f, map, z, arrays, (); workgroupsize)
end

#=
# Multi-GPU usage pattern (requires CUDA.jl):

# 1. Split the global mask and arrays into one chunk per device.
masks   = split_domain(global_mask, num_gpus)    # user-defined decomposition
X_devs  = split_domain(X, num_gpus)
Y_devs  = split_domain(Y, num_gpus)
Z_devs  = split_domain(Z, num_gpus)

# 2. Build one ActiveCellsMap per device, with indices living on that device.
maps = map(enumerate(masks)) do (i, m)
    CUDA.device!(i - 1)
    ActiveCellsMap(adapt(CUDABackend(), findall(m)))
end

# 3. Launch all devices concurrently; each apply! synchronizes its own stream.
@sync for (i, (acm, X_i, Y_i, Z_i)) in enumerate(zip(maps, X_devs, Y_devs, Z_devs))
    Threads.@spawn begin
        CUDA.device!(i - 1)
        apply!(f, acm, Z_i, (X_i, Y_i), p)
    end
end

# Notes:
#   - Domain decomposition and halo exchange across device boundaries are the
#     caller's responsibility; apply! itself is single-device.
#   - For stencil kernels, exchange ghost cells between devices before each apply! call.
=#
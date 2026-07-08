"""
$(TYPEDSIGNATURES)

Abstract supertype for the spatial discretizations of an `IceSheet`. Concrete subtypes such
as [`CommonGrid`](@ref) select, via multiple dispatch, how a component's arrays are laid out
and which coordinates they are evaluated on.
"""
abstract type AbstractGrid end

"""
$(TYPEDSIGNATURES)

Sentinel grid type: reuses the `common_grid` of the parent `IceSheet`,
avoiding redundant array allocation for components that share the same spatial discretization.
"""
struct CommonGrid <: AbstractGrid end

"""
$(TYPEDSIGNATURES)

A modular grid with uniform spacing and Cartesian or projected coordinate systems.
Supports 2D (nz = 1) and 3D configurations.

# Fields
- `nx`, `ny`, `nz`: number of cells in x, y, z directions.
- `x`, `y`, `z`: cell-centre coordinate vectors.
- `dx`, `dy`, `dz`: uniform cell spacings (scalars). Physical distances vary via `distortion`.
- `Lon`, `Lat`: geographic longitude and latitude (degrees) at each cell centre.
- `area`: horizontal cell area (m²), precomputed from `dx`, `dy`, and `distortion`.
- `distortion`: map-scale factor K; physical distance = `K * dx` (1 everywhere for Cartesian grids).
- `basins`: integer mask identifying drainage basins.
- `regions`: integer mask identifying user-defined regions.

Convenience constructors for a flat, regular, Cartesian 2D grid:

```julia
grid = RegularGrid(Float64, 6000e3, 6000e3, 16e3, 16e3)                       # CPU
grid = RegularGrid(CUDABackend(), Float32, 6000e3, 6000e3, 16e3, 16e3)        # GPU
```
"""
struct RegularGrid{I, T, V, M, MI} <: AbstractGrid
    nx::I
    ny::I
    nz::I
    x::V
    y::V
    z::V
    dx::T
    dy::T
    dz::T
    Lon::M
    Lat::M
    area::M
    distortion::M
    basins::MI
    regions::MI
end

function RegularGrid(T::Type{<:AbstractFloat}, lx, ly, dx, dy)
    return RegularGrid(CPU(), T, lx, ly, dx, dy)
end

function RegularGrid(backend::Backend, T::Type{<:AbstractFloat}, lx, ly, dx, dy)
    x_cpu = collect(range(T(0), step = T(dx), stop = T(lx))) .- T(lx) / 2
    y_cpu = collect(range(T(0), step = T(dy), stop = T(ly))) .- T(ly) / 2
    nx    = length(x_cpu)
    ny    = length(y_cpu)
    nz    = 1
    x     = KernelAbstractions.adapt(backend, x_cpu)
    y     = KernelAbstractions.adapt(backend, y_cpu)
    z     = KernelAbstractions.adapt(backend, [T(0)])
    Lon   = KernelAbstractions.zeros(backend, T, nx, ny)
    Lat   = KernelAbstractions.zeros(backend, T, nx, ny)
    area  = KernelAbstractions.adapt(backend, fill(T(dx * dy), nx, ny))
    dist  = KernelAbstractions.adapt(backend, ones(T, nx, ny))
    bas   = KernelAbstractions.zeros(backend, Int, nx, ny)
    reg   = KernelAbstractions.zeros(backend, Int, nx, ny)
    return RegularGrid(nx, ny, nz, x, y, z, T(dx), T(dy), T(1), Lon, Lat, area, dist, bas, reg)
end

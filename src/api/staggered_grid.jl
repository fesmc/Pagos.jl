"""
$(TYPEDSIGNATURES)

A [`Chmy.jl`](https://github.com/PTsolvers/Chmy.jl)-backed staggered (Arakawa C) grid.
Horizontal axes (`x`, `y`) are always uniform. The vertical axis is either:

 - a size-1 placeholder, for depth-integrated momentum balances
   (SIA, SSA). This mirrors [`RegularGrid`](@ref)'s `nz == 1` convention deliberately:
   [`MechanicState`](@ref)'s column (`M3`) fields rely on the grid always carrying a
   z-axis so that tensor/stress kernels stay dynamics-independent across the whole
   [`AbstractMomentumBalance`](@ref) hierarchy, with no 2D/3D special-casing.
 - a sigma-level `Chmy.FunctionAxis`, for full-column momentum balances
   (Blatter-Pattyn, Stokes), built from a [`CorrectedVerticalLayering`](@ref).

!!! note "Use `CorrectedVerticalLayering`, not `VerticalLayering`"
    Chmy derives an axis's cell centres *from* its vertices (`center = (vertex[i] +
    vertex[i+1]) / 2`) — the same face-first convention as [`CorrectedVerticalLayering`](@ref),
    so `zcenters(grid) == layering.ζ_aa` bit-for-bit. [`VerticalLayering`](@ref) goes
    the other way (vertices derived from independently-chosen midpoints) and does
    *not* round-trip through a `FunctionAxis` — passing it here would silently place
    layer midpoints where the grid does not think they are.

`StaggeredGrid` is a new, additive type: it coexists with [`RegularGrid`](@ref) rather
than replacing it (see `roadmaps/chmy.md`, Phase 1). `nx`, `ny`, `nz`, `x`, `y`, `z`,
`dx`, `dy`, `dz` are available as properties, computed from the wrapped Chmy grid, so
existing destructuring patterns like `(; nx, ny, nz) = grid` keep working. `dz` is only
defined when `nz == 1` — a full column's vertical spacing is non-uniform by
construction (see `ζ_aa`/`ζ_ac` on the `CorrectedVerticalLayering` instead).

!!! warning "Not node-for-node compatible with `RegularGrid`"
    Both convenience constructors take the same `(lx, ly, dx, dy)` inputs and both
    centre the *domain* on the origin, but they discretize it differently.
    `RegularGrid` places `lx/dx + 1` nodes spanning `[-lx/2, lx/2]` inclusive;
    `StaggeredGrid` follows Chmy's C-grid convention of `lx/dx` *cells*, whose
    centres run from `-lx/2 + dx/2` to `lx/2 - dx/2` (vertices sit on `±lx/2`).
    Same inputs therefore give `nx` smaller by one and cell-centre coordinates
    shifted by `dx/2` — anything sized against a `RegularGrid` (arrays, reference
    data) must be re-dimensioned when migrating, not just re-typed.

Unlike `RegularGrid`, `x`/`y`/`z` are lazy `Chmy` coordinate ranges (`xcenters`,
`ycenters`, `zcenters`), not arrays adapted onto `arch`'s device: Chmy kernels take a
grid and per-point indices, not a materialized coordinate array, so there is nothing
to adapt. `x`/`y`/`z` remain useful for host-side inspection, I/O and plotting.

# Fields
 - `arch`: the Chmy `Architecture` (device + execution backend) the grid was built for.
 - `grid`: the underlying Chmy `StructuredGrid{3}`.
 - `grid2d`: the same horizontal axes with a size-1 vertical axis — the grid the
   depth-integrated / depth-averaged fields of the state structs are dimensioned on
   (see below). Identical (`===`) to `grid` when `nz == 1`.
 - `Lon`, `Lat`: geographic longitude/latitude (degrees) at each horizontal cell centre.
 - `area`: horizontal cell area (m²), precomputed from `dx`, `dy`, and `distortion`.
 - `distortion`: map-scale factor K; physical distance = `K * dx` (1 everywhere for Cartesian grids).
 - `basins`: integer mask identifying drainage basins.
 - `regions`: integer mask identifying user-defined regions.

# Convenience constructors

```julia
grid = StaggeredGrid(Float64, 6000e3, 6000e3, 16e3, 16e3)                       # CPU, depth-integrated
grid = StaggeredGrid(CUDABackend(), Float32, 6000e3, 6000e3, 16e3, 16e3)        # GPU, depth-integrated
grid = StaggeredGrid(Float64, 6000e3, 6000e3, 16e3, 16e3, layering)            # full column
```

`topology` (a 2-tuple of `Chmy.Connectivity`, one per horizontal axis, applied to both
sides) defaults to `(Bounded(), Bounded())`. The z-axis connectivity is not user-settable:
it is always `Bounded()`, including for the size-1 depth-integrated placeholder.

!!! warning "Only `Bounded` and `Connected` are usable connectivities"
    Chmy 0.1.26 declares and exports four `Connectivity` types, but defines methods for
    only two of them. `Periodic` and `Flat` have *no methods anywhere in the package*, so
    rather than short-circuiting anything they fall through every dispatch:
    `Chmy.BoundaryConditions.batch_impl` covers `Bounded` and `Connected` only, and `bc!`
    batches over *all* axes, so a single `Periodic` or `Flat` axis turns `bc!` into a
    `MethodError` for every field on the grid — the `Bounded` axes included, and the
    depth-integrated SIA/SSA case along with them. Both are therefore rejected at
    construction rather than at the first boundary condition; periodic domains need
    periodic halo filling implemented Pagos-side first. This is also why the
    depth-integrated placeholder z-axis is `Bounded` and not the `Flat` that would
    describe it: `Bounded` costs a pair of never-read z ghost cells and keeps `bc!`
    working.

!!! note "Why there is a second grid (`grid2d`)"
    A Chmy `Field` takes its size from the grid it is built on (`size(grid, loc)`), so on
    a column grid *every* field is `nz`-deep. The state structs mix column fields
    (velocity, viscosity, the strain-rate and stress tensors) with genuinely
    depth-integrated ones (driving stress, depth-averaged velocity, basal friction, all
    of the topography), and the latter cannot be expressed as a location on the column
    grid: a `Vertex()` z-location gives `nz + 1` layers, not one. They are therefore
    built on `grid2d`, a second Chmy grid sharing the horizontal axes and carrying a
    size-1 `Bounded` z-axis. It is still a 3-dimensional grid — its fields are indexed
    `f[i, j, 1]`, never `f[i, j]`.
"""
struct StaggeredGrid{A,G,G2,M,MI} <: AbstractGrid
    arch::A
    grid::G
    grid2d::G2
    Lon::M
    Lat::M
    area::M
    distortion::M
    basins::MI
    regions::MI
end

const _StaggeredGridProperties = (:nx, :ny, :nz, :x, :y, :z, :dx, :dy, :dz)

function Base.getproperty(g::StaggeredGrid, s::Symbol)
    if s in (:arch, :grid, :grid2d, :Lon, :Lat, :area, :distortion, :basins, :regions)
        return getfield(g, s)
    end
    grid = getfield(g, :grid)
    s === :nx && return size(grid, Center())[1]
    s === :ny && return size(grid, Center())[2]
    s === :nz && return size(grid, Center())[3]
    s === :x && return xcenters(grid)
    s === :y && return ycenters(grid)
    s === :z && return zcenters(grid)
    s === :dx && return Δx(grid, Center(), 1, 1, 1)
    s === :dy && return Δy(grid, Center(), 1, 1, 1)
    if s === :dz
        size(grid, Center())[3] == 1 || error(
            "`dz` is not a single scalar for a column StaggeredGrid (the sigma " *
            "z-axis is non-uniform) — use the `ζ_aa`/`ζ_ac` vectors of the " *
            "`CorrectedVerticalLayering` it was built from instead.",
        )
        return Δz(grid, Center(), 1, 1, 1)
    end
    error("StaggeredGrid has no property $s")
end

Base.propertynames(::StaggeredGrid) = (
    :arch,
    :grid,
    :grid2d,
    :Lon,
    :Lat,
    :area,
    :distortion,
    :basins,
    :regions,
    _StaggeredGridProperties...,
)

# A `Connectivity` that Chmy declares but never dispatches on is not inert — it falls
# through dispatch instead of short-circuiting it. `BoundaryConditions.batch_impl` has
# methods for `Bounded` and `Connected` only, and `bc!` batches over *every* axis, so a
# single `Periodic` (or `Flat`) axis turns `bc!` into a `MethodError` for every field on
# the grid, the `Bounded` axes included. Rejecting it here trades a silent failure at the
# first boundary condition for a loud one at construction.
_check_connectivity(::Bounded, ::AbstractString) = nothing
_check_connectivity(::Connected, ::AbstractString) = nothing
_check_connectivity(conn::Connectivity, dim::AbstractString) = error(
    "`$(nameof(typeof(conn)))()` is not a usable topology for the $dim-axis of a " *
    "StaggeredGrid. Chmy $(pkgversion(Chmy)) declares and exports it but defines " *
    "no methods for it — in particular `BoundaryConditions.batch_impl` covers only " *
    "`Bounded` and `Connected`, and `bc!` batches over every axis, so one such axis " *
    "makes `bc!` a MethodError for every field on the grid (the `Bounded` axes " *
    "included). Use `Bounded()`, or `Connected()` for an MPI-exchanged axis.",
)

function _expand_topology(topology::NTuple{2,Connectivity})
    _check_connectivity(topology[1], "x")
    _check_connectivity(topology[2], "y")
    return (topology[1], topology[1]), (topology[2], topology[2])
end

"""
$(TYPEDSIGNATURES)

Build the sigma-level `Chmy.FunctionAxis` of a full-column [`StaggeredGrid`](@ref) from
the interface positions `layering.ζ_ac`.

The obvious one-liner — `FunctionAxis(i -> T(layering.ζ_ac[i]), layering.n)` — is wrong
twice over, and both failures are silent:

 - **Chmy evaluates the vertex function outside `1:n+1`.** `spacing(ax, Vertex(), i)` is
   `center(ax, i) - center(ax, i - 1)`, so `Δz`/`∂z` at the *bed* interface (`k == 1`)
   reaches `vertex(ax, 0)` and at the *surface* interface (`k == nz + 1`) reaches
   `vertex(ax, nz + 2)`. Those are interior points of every z-`Vertex` field — `ε̇_xz`,
   `ε̇_yz`, `w` — i.e. exactly the quantities a DIVA/Blatter-Pattyn balance is solved for.
   Chmy's `Launcher` additionally sweeps one halo ring (`Offset(-1)`), so even a
   cell-centred kernel reaches `i == 0`. Indexing `ζ_ac` directly is therefore an
   out-of-range read: a `BoundsError` in the best case, and under the `@inbounds` that
   kernels run with, silent garbage. The vertex function is made *total* instead,
   extrapolating linearly from the end spacings — the ghost interfaces this invents are
   never physical, but they make `Δz` at the bed and the surface well-defined.
 - **A closure over `layering` is not `isbits`**, because `CorrectedVerticalLayering`
   holds `Vector`s, and Chmy defines no `Adapt` rule for `StructuredGrid`/`FunctionAxis`
   to repair that downstream — so such a grid could never be passed to a GPU kernel. The
   interface positions are captured as an `NTuple` instead, which keeps the axis, and
   hence the whole grid, `isbits`.

The `NTuple` is passed by value into every kernel that takes the grid, so a very deep
column costs kernel-argument space (`8 * (nz + 1)` bytes); at the tens of layers ice-sheet
models use, that is negligible.

Cell centres are unaffected: Chmy derives them as `(ζ_ac[i] + ζ_ac[i+1]) / 2`, which
reproduces `layering.ζ_aa` bit-for-bit in both `Float32` and `Float64`.
"""
function _sigma_axis(T, layering::CorrectedVerticalLayering)
    n = layering.n
    ζ = ntuple(i -> T(layering.ζ_ac[i]), n + 1)
    Δlo = ζ[2] - ζ[1]
    Δhi = ζ[n+1] - ζ[n]
    # Branch-free: every arm is evaluated, so the index is always in range.
    zvertex(i) = ifelse(
        i < 1,
        ζ[1] + (i - 1) * Δlo,
        ifelse(i > n + 1, ζ[n+1] + (i - n - 1) * Δhi, ζ[clamp(i, 1, n + 1)]),
    )
    return FunctionAxis(zvertex, n)
end

function _geo_metadata(arch, T, nx, ny)
    backend = get_backend(arch)
    Lon = KernelAbstractions.zeros(backend, T, nx, ny)
    Lat = KernelAbstractions.zeros(backend, T, nx, ny)
    dist = KernelAbstractions.ones(backend, T, nx, ny)
    bas = KernelAbstractions.zeros(backend, Int, nx, ny)
    reg = KernelAbstractions.zeros(backend, Int, nx, ny)
    return Lon, Lat, dist, bas, reg
end

"""
$(TYPEDSIGNATURES)

Build a depth-integrated (`nz == 1`) `StaggeredGrid` on `arch`, spanning `(lx, ly)` at
uniform spacing `(dx, dy)`, centred on the origin. `grid2d === grid` here: there is only
one vertical layer, so depth-integrated and column fields live on the same grid.
"""
function StaggeredGrid(
    arch::Architecture,
    T::Type{<:AbstractFloat},
    lx,
    ly,
    dx,
    dy;
    topology::NTuple{2,Connectivity} = (Bounded(), Bounded()),
)
    nx = round(Int, lx / dx)
    ny = round(Int, ly / dy)
    xytopo = _expand_topology(topology)
    # `Bounded`, not the `Flat` that would describe this axis — see `_check_connectivity`.
    ztopo = (Bounded(), Bounded())

    xax = UniformAxis(T(-lx / 2), T(lx), nx)
    yax = UniformAxis(T(-ly / 2), T(ly), ny)
    zax = UniformAxis(T(0), T(1), 1)

    grid = StructuredGrid{typeof((xytopo..., ztopo))}(arch, xax, yax, zax)
    Lon, Lat, dist, bas, reg = _geo_metadata(arch, T, nx, ny)
    area = KernelAbstractions.adapt(get_backend(arch), fill(T(dx * dy), nx, ny))
    return StaggeredGrid(arch, grid, grid, Lon, Lat, area, dist, bas, reg)
end

StaggeredGrid(T::Type{<:AbstractFloat}, lx, ly, dx, dy; kwargs...) =
    StaggeredGrid(Arch(KernelAbstractions.CPU()), T, lx, ly, dx, dy; kwargs...)

StaggeredGrid(
    backend::KernelAbstractions.Backend,
    T::Type{<:AbstractFloat},
    lx,
    ly,
    dx,
    dy;
    kwargs...,
) = StaggeredGrid(Arch(backend), T, lx, ly, dx, dy; kwargs...)

"""
$(TYPEDSIGNATURES)

Build a full-column `StaggeredGrid` on `arch`: uniform `(lx, ly)`/`(dx, dy)` horizontal
axes plus a sigma-level vertical axis (`Chmy.FunctionAxis`) built from `layering`
(see the note on [`CorrectedVerticalLayering`](@ref) vs [`VerticalLayering`](@ref)
above — only the former round-trips through Chmy's axis).
"""
function StaggeredGrid(
    arch::Architecture,
    T::Type{<:AbstractFloat},
    lx,
    ly,
    dx,
    dy,
    layering::CorrectedVerticalLayering;
    topology::NTuple{2,Connectivity} = (Bounded(), Bounded()),
)
    nx = round(Int, lx / dx)
    ny = round(Int, ly / dy)
    xytopo = _expand_topology(topology)
    ztopo = (Bounded(), Bounded())

    xax = UniformAxis(T(-lx / 2), T(lx), nx)
    yax = UniformAxis(T(-ly / 2), T(ly), ny)
    zax = _sigma_axis(T, layering)

    topo = typeof((xytopo..., ztopo))
    grid = StructuredGrid{topo}(arch, xax, yax, zax)
    # Same horizontal axes, single vertical layer: the grid that carries the
    # depth-integrated fields of the state structs (see the `grid2d` note above).
    grid2d = StructuredGrid{topo}(arch, xax, yax, UniformAxis(T(0), T(1), 1))
    Lon, Lat, dist, bas, reg = _geo_metadata(arch, T, nx, ny)
    area = KernelAbstractions.adapt(get_backend(arch), fill(T(dx * dy), nx, ny))
    return StaggeredGrid(arch, grid, grid2d, Lon, Lat, area, dist, bas, reg)
end

StaggeredGrid(
    T::Type{<:AbstractFloat},
    lx,
    ly,
    dx,
    dy,
    layering::CorrectedVerticalLayering;
    kwargs...,
) = StaggeredGrid(
    Arch(KernelAbstractions.CPU()),
    T,
    lx,
    ly,
    dx,
    dy,
    layering;
    kwargs...,
)

StaggeredGrid(
    backend::KernelAbstractions.Backend,
    T::Type{<:AbstractFloat},
    lx,
    ly,
    dx,
    dy,
    layering::CorrectedVerticalLayering;
    kwargs...,
) = StaggeredGrid(Arch(backend), T, lx, ly, dx, dy, layering; kwargs...)

"""
$(TYPEDSIGNATURES)

Fill `f` from an analytic function of coordinates, `fun(x, y, ζ) -> value` (or, with
`discrete = true`, `fun(grid, loc, i, j, k) -> value`), against the `StaggeredGrid` `f`
was built on — extending `Chmy.set!` so callers never need `f`'s bare `.grid`/`.grid2d`.

Dispatches to whichever of `sg.grid`/`sg.grid2d` actually matches `f`'s own z-extent
(`size(f, 3) == 1` selects `sg.grid2d`), rather than trusting the caller to pick the
right one. That trust would be misplaced: both are 3D `Chmy.StructuredGrid`s (`grid2d`
carries a size-1 z-axis, not no z-axis, so its own `Field`s are also 3D — see the
`StaggeredGrid` docstring), so passing the wrong one is not a `MethodError` — it is a
silently wrong third coordinate, since each grid's z-axis answers `coord` independently.

!!! warning "The third coordinate is the sigma level ζ, not physical elevation z"
    On a full-column grid, `sg.grid`'s z-axis is [`CorrectedVerticalLayering`](@ref)'s
    dimensionless sigma level (see `_sigma_axis`), so `fun`'s third argument is `ζ ∈
    [0, 1]`, not a physical depth. An analytic function written against physical depth
    must convert inside `fun` using the local ice thickness — nothing here does that
    conversion.

# Examples

```jldoctest
julia> sg = StaggeredGrid(Float64, 4.0, 4.0, 1.0, 1.0);

julia> f = Field(sg.arch, sg.grid, Center());

julia> set!(f, sg, (x, y, ζ) -> x + y);

julia> interior(f)[1, 1, 1] ≈ sg.x[1] + sg.y[1]
true
```
"""
function Chmy.set!(f::Field{T,3}, sg::StaggeredGrid, fun; kwargs...) where {T}
    grid = size(interior(f), 3) == 1 ? sg.grid2d : sg.grid
    return Chmy.set!(f, grid, fun; kwargs...)
end

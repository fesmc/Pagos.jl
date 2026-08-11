"""
$(TYPEDSIGNATURES)

The execution context a Chmy-native Pagos kernel needs: *where* to run (`arch`), *over
what* (`grid`/`grid2d`), and *how to launch* (`launch`/`launch2d`). One `Runtime` is built
per [`StaggeredGrid`](@ref) and threaded through the `Field`-based dispatch of the physics
signatures (`creep!(cf::AbstractField, σ_e, law, rt)`, ...).

The plain-array dispatch of those same functions (`creep!(cf::AbstractArray, σ_e, law)`)
takes **no** `Runtime`: it reads its backend off the array itself
(`KernelAbstractions.get_backend`) and has no grid, halo, or launcher concept.

# Fields
 - `arch`: the Chmy `Architecture` (device + backend); taken from the grid.
 - `grid`: the **bare Chmy grid** (`StructuredGrid{3}`), not the Pagos
   [`StaggeredGrid`](@ref) it was built from — Chmy's grid operators (`∂x`, `Δx`, ...),
   `Launcher` and `bc!` all dispatch on `Chmy.StructuredGrid`, and the wrapper has no
   `Adapt` rule, so it can never be passed into a kernel. Reach the geographic metadata
   (`area`, `distortion`, ...) through [`IceSheet`](@ref)'s [`StaggeredGrid`](@ref) and
   pass the array into the kernel explicitly.
 - `grid2d`: the bare Chmy grid the *depth-integrated* fields live on (`StaggeredGrid`'s
   `grid2d`: same horizontal axes, size-1 z-axis). `=== grid` when `nz == 1`.
 - `launch`: a Chmy `Launcher` sized for `grid`. Sweeps `size(grid, Center()) .+ 2` with
   an `Offset(-1)`, i.e. one halo ring beyond the interior, so a kernel launched this way
   also fills that ring.
 - `launch2d`: a [`FlatLauncher`](@ref) sized for `grid2d` — same call signature as
   `Launcher`, minus the phantom z ring Chmy's worksize rule adds on a size-1 axis. Falls
   back to a plain `Launcher` when `outer_width` is set.

# Launching

`Runtime` adds no launching API of its own; the fields are exactly the three arguments
Chmy's `Launcher` wants, so launch the way Chmy documents it:

```julia
rt = Runtime(grid)
rt.launch(rt.arch, rt.grid, my_kernel! => (out, in, rt.grid))
rt.launch(rt.arch, rt.grid, my_kernel! => (out, in, rt.grid); bc = batch(rt.grid, out => Neumann()))
```

!!! warning "Pick the launcher that matches the *output* field's grid"
    A `Field`'s size comes from the grid it was built on, and the state structs mix both
    (see [`StaggeredGrid`](@ref)'s `grid2d` note). A mismatch does **not** error: it
    silently sweeps too few layers (2D launcher over a column field) or runs off the end
    of the shallow field's `k` range (column launcher over a depth-integrated one). A
    kernel that reads both — SIA/SSA driving stress from a column viscosity, DIVA's
    vertical integrals — launches on whichever grid its *output* lives on and indexes the
    other explicitly.

!!! warning "`outer_width` is unsafe under AD"
    It enables Chmy's communication/computation overlap: the launcher splits the sweep
    into an interior part and boundary slabs run on async `Worker` Tasks. Keep
    `outer_width = nothing` inside an Enzyme-differentiated region.
"""
struct Runtime{A,G,G2,L,L2}
    arch::A
    grid::G
    grid2d::G2
    launch::L
    launch2d::L2
end

"""
$(TYPEDSIGNATURES)

A drop-in `Chmy.Launcher` for **depth-integrated** grids: same call signature, same
`Offset` convention in `x`/`y`, but worksize `(nx + 2, ny + 2, 1)` with `Offset(-1, -1, 0)`
instead of Chmy's `size(grid, Center()) .+ 2` = `(nx + 2, ny + 2, 3)` with `Offset(-1)`.

`grid2d`'s z axis has extent 1, so the two extra planes Chmy's rule adds are a halo ring in
a dimension that has none, and nothing can observe them: the 2D grid operators never index
`k ± 1`, column kernels read a depth-integrated field at an explicit `k = 1`, and `interior`
on a `grid2d` `Field` returns `k = 1` only. Skipping them measures **2.22× on CPU** and
**1.99× on GPU** (`benchmark/basics/gpu/README.md` §1).

Worksize *and* groupsize are type parameters rather than fields, which is what lets a
static-size KA kernel specialise its index arithmetic and bounds checks at compile time.

!!! warning "Depth-integrated grids only"
    Construction throws unless `size(grid, Center())[3] == 1`. Column kernels keep the
    ordinary `Launcher` through `rt.launch`, which still sweeps its z ring — a real one.

`outer_width` is not implemented here; [`Runtime`](@ref) falls back to a plain `Launcher`
for `launch2d` when it is set. The `synchronize` after each launch is kept, matching
`Launcher` (`benchmark/basics/gpu/README.md` §2).
"""
struct FlatLauncher{Worksize,GroupSize,B}
    backend::B
end

function FlatLauncher(arch, grid)
    n = size(grid, Center())
    n[3] == 1 || throw(
        ArgumentError(
            "FlatLauncher needs a depth-integrated grid (z extent 1), got nz = $(n[3]). " *
            "Column kernels launch with `rt.launch`, not `rt.launch2d`.",
        ),
    )
    backend = get_backend(arch)
    groupsize = heuristic_groupsize(backend, Val(3))
    ws = (n[1] + 2, n[2] + 2, 1)
    return FlatLauncher{ws,groupsize,typeof(backend)}(backend)
end

# Extend Chmy's own generics, not new same-named ones: `worksize`/`outer_width` reach Pagos
# only through `@reexport using Chmy`, so a bare definition here would shadow Chmy's and
# make the unqualified name ambiguous at every call site.
Base.@assume_effects :foldable Base.ndims(::FlatLauncher{WS}) where {WS} = length(WS)
Base.@assume_effects :foldable Chmy.KernelLaunch.worksize(::FlatLauncher{WS}) where {WS} =
    WS
Base.@assume_effects :foldable Chmy.KernelLaunch.outer_width(::FlatLauncher) = nothing

# `GS`/`WS` are passed to `kernel(backend, groupsize, ndrange)` from the type, not from a
# field: a field is read at runtime, so KA's `StaticSize` wrapper cannot fold it and the
# launch config never reaches the type system.
function (launcher::FlatLauncher{WS,GS})(
    arch::Architecture,
    grid,
    kernel_and_args::Pair{F,Args};
    bc = nothing,
) where {WS,GS,F,Args}
    kernel, args = kernel_and_args
    offset = Offset(-1, -1, 0)

    if isnothing(bc)
        kernel(launcher.backend, GS, WS)(args..., offset)
    else
        # Mirrors Chmy's own `launch_with_bc` on the `outer_width === nothing` branch:
        # whole-domain kernel first, then one `bc!` over the batch.
        gs = KernelAbstractions.NDIteration.StaticSize(GS)
        ws = KernelAbstractions.NDIteration.StaticSize(WS)
        kernel(launcher.backend, gs, ws)(args..., offset)
        bc!(arch, grid, bc)
    end

    KernelAbstractions.synchronize(launcher.backend)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Build a [`Runtime`](@ref) for `grid`, reusing the architecture the grid was constructed
on. `outer_width` is forwarded to Chmy's `Launcher` (see the [`Runtime`](@ref) docstring
on when *not* to set it).

`grid` gets an ordinary Chmy `Launcher`; `grid2d` gets a [`FlatLauncher`](@ref), even when
`nz == 1` and the two grids are the *same object* — a `FlatLauncher` owns no `Worker` Tasks,
so building a second one costs nothing. With `outer_width` set, `launch2d` falls back to a
plain `Launcher`, and *there* the same-object case reuses `launch` rather than spawning a
duplicate set of `Worker` Tasks.
"""
function Runtime(grid::StaggeredGrid; outer_width = nothing)
    arch = grid.arch
    g3, g2 = grid.grid, grid.grid2d
    launch = Launcher(arch, g3; outer_width)
    launch2d = if isnothing(outer_width)
        FlatLauncher(arch, g2)
    else
        g2 === g3 ? launch : Launcher(arch, g2; outer_width)
    end
    return Runtime(arch, g3, g2, launch, launch2d)
end

KernelAbstractions.get_backend(rt::Runtime) = KernelAbstractions.get_backend(rt.arch)

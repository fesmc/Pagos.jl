"""
$(TYPEDSIGNATURES)

The execution context a Chmy-native Pagos kernel needs: *where* to run (`arch`), *over
what* (`grid`/`grid2d`), and *how to launch* (`launch`/`launch2d`). One `Runtime` is built
per [`StaggeredGrid`](@ref) and threaded through the `Field`-based dispatch of the physics
signatures (`creep!(cf::AbstractField, σ_e, law, rt)`, ...).

The plain-array dispatch of those same functions (`creep!(cf::AbstractArray, σ_e, law)`)
takes **no** `Runtime`: it reads its backend off the array itself
(`KernelAbstractions.get_backend`) and has no grid, halo, or launcher concept. `Runtime`
is internal machinery — it never appears in a signature a user is expected to call with
plain arrays (see `pagos-roadmaps/chmy.md`, Phase 2).

# Fields
 - `arch`: the Chmy `Architecture` (device + backend); taken from the grid.
 - `grid`: the **bare Chmy grid** (`StructuredGrid{3}`), not the Pagos
   [`StaggeredGrid`](@ref) it was built from — see below.
 - `grid2d`: the bare Chmy grid the *depth-integrated* fields live on (`StaggeredGrid`'s
   `grid2d`: same horizontal axes, size-1 z-axis). `=== grid` when `nz == 1`.
 - `launch`: a Chmy `Launcher` sized for `grid`.
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
    (see [`StaggeredGrid`](@ref)'s `grid2d` note). Launch a kernel writing
    depth-integrated fields (ice thickness, driving stress, mass fluxes) with
    `rt.launch2d(rt.arch, rt.grid2d, ...)`, and one writing column fields with
    `rt.launch(rt.arch, rt.grid, ...)`. A mismatch does **not** error: it silently sweeps
    too few layers (2D launcher over a column field) or runs off the end of the shallow
    field's `k` range (column launcher over a depth-integrated one). A kernel that reads
    both — SIA/SSA driving stress from a column viscosity, DIVA's vertical integrals —
    launches on whichever grid its *output* lives on and indexes the other explicitly.

!!! note "`rt.grid` is the Chmy grid, not the `StaggeredGrid`"
    Chmy's grid operators (`∂x`, `Δx`, ...), `Launcher` and `bc!` all dispatch on
    `Chmy.StructuredGrid`; the Pagos [`StaggeredGrid`](@ref) wrapper is opaque to them
    and has no `Adapt` rule, so it can never be passed into a kernel. Holding the bare
    grid keeps `Runtime` usable directly with Chmy's documented call signature. The
    geographic metadata on [`StaggeredGrid`](@ref) (`area`, `distortion`, ...) is
    reached through the grid itself, which [`IceSheet`](@ref) owns; a kernel that needs
    it takes the array as an explicit argument, as Chmy kernels do for any other data.

!!! note "Worksize includes one halo ring"
    Chmy's `Launcher` sweeps `size(grid, Center()) .+ 2` points with an `Offset(-1)`,
    i.e. one halo ring beyond the interior, so a kernel launched this way also fills
    that ring. Whether that ring should be *computed* or `bc!`-filled is a per-kernel
    decision (see `pagos-roadmaps/chmy.md`, Phase 4).

    `launch2d` is the exception, and deliberately so: it is a [`FlatLauncher`](@ref),
    which sweeps the ring in `x`/`y` but **not** in `z`, because `grid2d`'s z axis has
    extent 1 and therefore has no ring to fill. See that docstring for why nothing reads
    the planes it stops writing.

# `outer_width` and AD

`outer_width` (default `nothing`) enables Chmy's communication/computation overlap:
the launcher splits the sweep into an interior part and boundary slabs run on async
`Worker` tasks. That path spawns Tasks and is **not** safe inside an
Enzyme-differentiated region — keep `outer_width = nothing` there (see
`pagos-roadmaps/chmy.md`, Phase 5).
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
a dimension that has none: every depth-integrated kernel reads and writes **three** k-planes
where only `k = 1` holds data. Removing them measures **2.22× on CPU** (4 threads, 380²,
Float64, 100 fixed PT iterations: 15.83 → 7.13 ms/iter, velocities bit-identical) and
**1.99× on GPU**, 2.9–4.2× per kernel (`benchmark/basics/gpu/README.md` §1).

# Why nothing reads the planes this stops writing

 - The 2D grid operators (`∂x`, `∂y`, `lerp`/`hlerp` at `aa`/`ab`/`acx`/`acy`) never index
   `k ± 1`, so no depth-integrated kernel can observe them.
 - Column kernels that read a depth-integrated field index it explicitly at `k = 1` — e.g.
   `_dotvel_staggered_bp!` writes `lerp(H, NODE_ACX, grid2d, i, j, 1)`.
 - `interior` on a `grid2d` `Field` returns `k = 1` only, so neither output nor tests see
   them.

Every kernel launched through `launch2d` takes `i, j` from the sweep and then either writes
a `grid2d` field at the sweep's own index (only `k = 1` meaningful) or writes a *column*
field through an explicit internal `for k` loop (`_depthaverage!`, `_velocities3D_ssa!`,
`_verticalvelocity!`, `_viscosity_integrals!`, `_vertical_line_relax_bp!`) — in which case
the sweep's `k` was pure repetition, running the same column integral three times.

!!! warning "Depth-integrated grids only"
    Construction throws unless `size(grid, Center())[3] == 1`. Column kernels keep the
    ordinary `Launcher` through `rt.launch`, which still sweeps its z ring — a real one.

!!! note "`outer_width` is not implemented here"
    The communication/computation overlap path (async `Worker` tasks, per-side `bc!`) is
    Chmy's; [`Runtime`](@ref) falls back to a plain `Launcher` for `launch2d` when
    `outer_width` is set, so that configuration is unchanged by this type.

The `synchronize` after each launch is kept, matching `Launcher`. Dropping it is worth a
further ~1.3× but is a Chmy-wide policy question and size-dependent (a loss at 381²) — see
`benchmark/basics/gpu/README.md` §2.
"""
struct FlatLauncher{Worksize,B,G}
    backend::B
    groupsize::G
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
    return FlatLauncher{ws,typeof(backend),typeof(groupsize)}(backend, groupsize)
end

# Extend Chmy's own generic functions, not new same-named ones: `worksize`/`outer_width`
# reach Pagos only through `@reexport using Chmy`, so a *bare* definition here would create
# a distinct `Pagos.worksize` that shadows Chmy's and makes the unqualified name ambiguous
# at every call site, including the existing `worksize(rt.launch)` ones.
Base.@assume_effects :foldable Base.ndims(::FlatLauncher{WS}) where {WS} = length(WS)
Base.@assume_effects :foldable Chmy.KernelLaunch.worksize(::FlatLauncher{WS}) where {WS} =
    WS
Base.@assume_effects :foldable Chmy.KernelLaunch.outer_width(::FlatLauncher) = nothing

function (launcher::FlatLauncher{WS})(
    arch::Architecture,
    grid,
    kernel_and_args::Pair{F,Args};
    bc = nothing,
) where {WS,F,Args}
    kernel, args = kernel_and_args
    offset = Offset(-1, -1, 0)

    if isnothing(bc)
        kernel(launcher.backend, launcher.groupsize, WS)(args..., offset)
    else
        # Mirrors Chmy's own `launch_with_bc` on the `outer_width === nothing` branch:
        # whole-domain kernel first, then one `bc!` over the batch.
        gs = KernelAbstractions.NDIteration.StaticSize(launcher.groupsize)
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

`grid` gets an ordinary Chmy `Launcher`; `grid2d` gets a [`FlatLauncher`](@ref), which
drops the phantom z ring Chmy's worksize rule adds on a size-1 axis. This holds even when
`nz == 1` and the two grids are the *same object*: the depth-integrated launcher is still
the one that should sweep one k-plane, and a `FlatLauncher` owns no `Worker` Tasks, so
building a second launcher costs nothing.

The exception is `outer_width`: that path is Chmy's async communication/computation overlap,
which [`FlatLauncher`](@ref) does not implement. With it set, `launch2d` falls back to a
plain `Launcher` — and *there* the same-object reuse still matters, since each `Launcher`
would otherwise spawn a duplicate set of `Worker` Tasks.
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

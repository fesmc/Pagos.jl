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
plain arrays (see `roadmaps/chmy.md`, Phase 2).

# Fields
 - `arch`: the Chmy `Architecture` (device + backend); taken from the grid.
 - `grid`: the **bare Chmy grid** (`StructuredGrid{3}`), not the Pagos
   [`StaggeredGrid`](@ref) it was built from — see below.
 - `grid2d`: the bare Chmy grid the *depth-integrated* fields live on (`StaggeredGrid`'s
   `grid2d`: same horizontal axes, size-1 z-axis). `=== grid` when `nz == 1`.
 - `launch`: a Chmy `Launcher` sized for `grid`.
 - `launch2d`: a Chmy `Launcher` sized for `grid2d` (`=== launch` when `nz == 1`).

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
    decision (see `roadmaps/chmy.md`, Phase 4).

# `outer_width` and AD

`outer_width` (default `nothing`) enables Chmy's communication/computation overlap:
the launcher splits the sweep into an interior part and boundary slabs run on async
`Worker` tasks. That path spawns Tasks and is **not** safe inside an
Enzyme-differentiated region — keep `outer_width = nothing` there (see
`roadmaps/chmy.md`, Phase 5).
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

Build a [`Runtime`](@ref) for `grid`, reusing the architecture the grid was constructed
on. `outer_width` is forwarded to Chmy's `Launcher` (see the [`Runtime`](@ref) docstring
on when *not* to set it).

Both of the grid's Chmy grids get a `Launcher`. When they are the same object (`nz == 1`,
where `StaggeredGrid.grid2d === .grid`) the *same* `Launcher` is reused rather than a
second one built: with `outer_width` set, a `Launcher` owns async `Worker` Tasks, so a
duplicate would spawn a second set of them for no benefit.
"""
function Runtime(grid::StaggeredGrid; outer_width = nothing)
    arch = grid.arch
    g3, g2 = grid.grid, grid.grid2d
    launch = Launcher(arch, g3; outer_width)
    launch2d = g2 === g3 ? launch : Launcher(arch, g2; outer_width)
    return Runtime(arch, g3, g2, launch, launch2d)
end

KernelAbstractions.get_backend(rt::Runtime) = KernelAbstractions.get_backend(rt.arch)

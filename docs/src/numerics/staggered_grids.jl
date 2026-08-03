#=

# [Staggered grids](@id staggered_grids)

Ice-dynamics solvers live on a **staggered (Arakawa C) grid**: scalars (thickness,
viscosity) sit at cell centres, while the two velocity components sit on the cell
edges, each offset by half a cell in one direction. Pagos builds this on top of
[Chmy.jl](https://github.com/PTsolvers/Chmy.jl), which Pagos re-exports in full —
every name below (`Center`, `Vertex`, `UniformGrid`, `Field`, `∂x`, `lerp`, ...) is
available directly after `using Pagos`, with no separate `using Chmy` needed.

Chmy's model rests on two ideas:

1. **Every `Field` carries a location** `(LX, LY[, LZ])`, each entry `Center()` or
   `Vertex()`. A vertex sits on the boundary between two centres, so a
   `Vertex`-located axis has one more point than its `Center` counterpart.
2. **Operators are explicit, index-level functions**, not a lazy expression tree:
   `∂x(f, grid, I...)` returns the value of the stencil at grid index `I`. There is
   no automatic interpolation when locations disagree — *you* choose the index at
   which to evaluate, and reconcile mismatched locations yourself with
   [`lerp`](@ref)/[`hlerp`](@ref) where the physics calls for it. This is the
   opposite trade-off from a lazy `@at`-style macro: less magic, but the location of
   every access is visible at the call site, and it composes cleanly with kernel
   compilation and, eventually, automatic differentiation.

## Locations and the Pagos node naming

The familiar Pagos / Yelmo node names map directly onto a pair of horizontal
locations:

| Pagos node | location `(LX, LY)`   | typical field                                     |
|:---------- |:---------------------- |:-------------------------------------------------- |
| `aa`       | `(Center, Center)`     | `H`, `z_srf`, `z_bed`, viscosity                   |
| `acx`      | `(Vertex, Center)`     | `u`, `taud_acx`, x-flux                            |
| `acy`      | `(Center, Vertex)`     | `v`, `taud_acy`, y-flux                            |
| `ab`       | `(Vertex, Vertex)`     | corner-staggered shear strain rate / stress        |

A [`VectorField`](@ref) built from a `Center` grid already places its `.x`/`.y`
components at `acx`/`acy` for you — this is the natural home for a velocity or a
flux.

## The third axis: layer midpoints vs layer interfaces

The vertical axis stages exactly the same way, in the terrain-following
``\sigma`` coordinate (``\sigma = 0`` at the bed, ``1`` at the surface):

- **layer midpoints** ``\zeta_{aa}`` are z-`Center` — where a layer *quantity* lives
  (viscosity, temperature, the horizontal velocity of a layer);
- **layer interfaces** ``\zeta_{ac}`` are z-`Vertex` — where anything *differentiated
  with respect to* ``z`` lives (the vertical shear ``\dot\varepsilon_{xz}``,
  ``\dot\varepsilon_{yz}``, and the vertical velocity ``w``).

A z-`Vertex` field therefore has `nz + 1` layers, not `nz`. In Pagos' node names the
suffix `_ac` marks the interface variants: `aa_ac`, `acx_ac`, `acy_ac`.

As in the horizontal, the choice is forced by the operators rather than free: `∂z` maps
`Center → Vertex`, so `∂u/∂z` with `u` at `acx` (z-`Center`) lands on `acx_ac`, which is
exactly where ``\dot\varepsilon_{xz} = (\partial u/\partial z + \partial w/\partial x)/2``
needs it — and `∂w/∂x` with `w` at `aa_ac` lands on the same `acx_ac`, so the two terms
add with no interpolation. The whole tensor closes this way, which is the check that the
layout is right.

!!! warning "On a non-uniform ``\sigma`` axis, use [`∂z_σ`](@ref), not Chmy's `∂z`"
    Chmy scales the vertical difference by the spacing at the *field's* location rather
    than the result's. On a `UniformAxis` the two coincide, so nothing in the horizontal
    is affected; on the stretched ``\sigma`` axis they differ, and `∂z` is wrong by tens
    of percent mid-column.

## [The `2`/`3` suffix: which grid a field lives on](@id node_type_parameters)

Pagos' state structs (`MechanicState`, `StressState`, ...) carry **one type parameter per
distinct node class**, so that a mis-wired constructor is a `MethodError` at construction
rather than a half-cell shift at runtime. Those parameter names are what you see in the
struct definitions, and they combine the location above with *which of the two grids* the
field is built on:

- a [`StaggeredGrid`](@ref) carries `grid` (the full column, `nz` layers) **and**
  `grid2d` (the same horizontal axes with a size-1 vertical axis);
- the trailing **`3`** means the column grid, the trailing **`2`** means `grid2d`, i.e.
  a depth-integrated quantity;
- a `Z` before the digit marks a z-`Vertex` (interface) field.

Concretely, for a grid with `nx = ny = 8`, `nz = 6` — every shape below is the actual
`size(interior(f))`:

| Parameter | Location `(LX, LY, LZ)` | Grid     | Shape       | Example field                 |
|:--------- |:----------------------- |:-------- |:----------- |:----------------------------- |
| `AA2`     | `(C, C, C)`             | `grid2d` | `(8, 8, 1)` | `topography.thickness`        |
| `ACX2`    | `(V, C, C)`             | `grid2d` | `(9, 8, 1)` | `velocity.depthaverage_x` (ū) |
| `ACY2`    | `(C, V, C)`             | `grid2d` | `(8, 9, 1)` | `velocity.depthaverage_y` (v̄) |
| `AB2`     | `(V, V, C)`             | `grid2d` | `(9, 9, 1)` | `velocity.depthaverage_x_dy`  |
| `AA3`     | `(C, C, C)`             | `grid`   | `(8, 8, 6)` | `material.viscosity` µ(z)     |
| `ACX3`    | `(V, C, C)`             | `grid`   | `(9, 8, 6)` | `velocity.x` — u(z)           |
| `ACY3`    | `(C, V, C)`             | `grid`   | `(8, 9, 6)` | `velocity.y` — v(z)           |
| `AB3`     | `(V, V, C)`             | `grid`   | `(9, 9, 6)` | `strainrate.xy`               |
| `AAZ3`    | `(C, C, V)`             | `grid`   | `(8, 8, 7)` | `velocity.z` — w              |
| `ACXZ3`   | `(V, C, V)`             | `grid`   | `(9, 8, 7)` | `strainrate.xz`               |
| `ACYZ3`   | `(C, V, V)`             | `grid`   | `(8, 9, 7)` | `strainrate.yz`               |

Reading a parameter name is therefore mechanical: `ACXZ3` = x-`Vertex`, y-`Center`,
z-`Vertex`, on the column grid — an x-face, layer-interface field.

Two consequences worth internalising, because both have caused real bugs:

1. **`AA2` and `AA3` are the same shape when `nz == 1`, and only then.** A depth-integrated
   grid has `grid2d === grid`, so the two collapse and code can read a column field as
   though it were 2D without complaint. That is why a solver written for `nz == 1` can
   silently depend on the collapse — and why it breaks the moment a real column appears.
2. **A z-`Vertex` field on a depth-integrated grid has *2* layers, not 1** — the bed and
   the surface. So `strainrate.xz` is *not* shaped like `strainrate.xx` even when `nz == 1`,
   and there is deliberately no `*Z2` parameter.

=#

using Pagos
using KernelAbstractions: @kernel, @index, CPU

arch = Arch(CPU())
grid = UniformGrid(arch; origin = (0.0, 0.0), extent = (1.0, 1.0), dims = (8, 8))
launch = Launcher(arch, grid)

## A scalar field on aa-nodes (cell centres).
H = Field(arch, grid, Center())
set!(H, grid, (x, y) -> x^2)
interior(H)

#=

## Staggering an aa-field onto the ac-nodes

A scalar known at cell centres (`aa`) staggered onto a velocity point (`acx`/`acy`)
is a half-cell average — [`lerp`](@ref). Unlike a lazy macro, this is an ordinary
kernel: you write the loop body once and [`Launcher`](@ref) dispatches it across
threads (CPU) or blocks (GPU).

=#

@kernel inbounds = true function stagger_kernel!(out, f, g, O)
    I = @index(Global, NTuple)
    I = I + O
    out[I...] = lerp(f, location(out), g, I...)
end

H_acx = Field(arch, grid, (Vertex(), Center()))
launch(arch, grid, stagger_kernel! => (H_acx, H, grid))
interior(H_acx)

#=

## Derivatives: aa → ac

The `x`-derivative of an `aa`-field naturally lives on the `acx` velocity points —
exactly where a flux-divergence solver wants it. This is the same idiom Chmy's own
test suite uses to build a gradient into a [`VectorField`](@ref):

=#

@kernel inbounds = true function grad_kernel!(q, H, g, O)
    I = @index(Global, NTuple)
    I = I + O
    q.x[I...] = ∂x(H, g, I...)
    q.y[I...] = ∂y(H, g, I...)
end

q = VectorField(arch, grid)
launch(arch, grid, grad_kernel! => (q, H, grid))

## `H = x²` so `∂H/∂x = 2x`; central differences are exact for a quadratic.
## The outermost ring of vertices/centres borders the domain edge, where `H`'s
## halo has not been filled by a boundary condition yet (see the note below) —
## compare the interior only:
interior(q.x)[2:end-1, 4]

#=

## A compound expression, evaluated back at the centres

`lerp` and `∂x`/`∂y` are ordinary functions, so composing them is just more code in
the same kernel — e.g. an effective-strain-rate-like quantity `(∂x H)² + (∂y H)²`,
staggered back onto the `aa` nodes it started from:

=#

@kernel inbounds = true function grad_squared_kernel!(out, q, g, O)
    I = @index(Global, NTuple)
    I = I + O
    qx_c = lerp(q.x, location(out), g, I...)
    qy_c = lerp(q.y, location(out), g, I...)
    out[I...] = qx_c^2 + qy_c^2
end

grad2 = Field(arch, grid, Center())
launch(arch, grid, grad_squared_kernel! => (grad2, q, grid))
interior(grad2)[2:end-1, 4]

#=

!!! note "Boundaries need an explicit boundary condition"
    Unlike the previous index-wrapping implementation, Chmy fields carry an
    explicit halo (`Field(...; halo=1)` by default) rather than wrapping indices at
    the domain edge. Reading from that halo before it has been filled — via
    [`bc!`](@ref) or a halo exchange — gives whatever was last written there, not a
    periodic/reflective/one-sided value; that is why the examples above only
    inspect interior points. Mapping Pagos' current boundary conventions
    (`FlatIndexing`, `ReflectiveIndexing`, `PeriodicIndexing`) onto `bc!` is tracked
    as its own step in `roadmaps/chmy.md`.

## Backends and distribution

Everything above ran on `CPU()`. The same kernels run unchanged on `CUDABackend()`,
`ROCBackend()`, or `MetalBackend()` by constructing `arch = Arch(backend)`
accordingly, and Chmy's `Distributed` module (`CartesianTopology`,
`exchange_halo!`) extends the same `Field`/`Launcher` model across MPI ranks for
multi-GPU runs — see the [Chmy.jl documentation](https://ptsolvers.github.io/Chmy.jl/)
for details.

## Summary

- `aa = (Center, Center)`, `acx = (Vertex, Center)`, `acy = (Center, Vertex)`,
  `ab = (Vertex, Vertex)`; vertically, layer midpoints are z-`Center` and layer
  interfaces z-`Vertex` (the `_ac` suffix), the latter with `nz + 1` layers.
- **Reading a state-struct type parameter** (`AA2`, `ACXZ3`, ...): the letters give the
  location, a `Z` marks z-`Vertex`, and the digit says which grid — `2` for the
  depth-integrated `grid2d`, `3` for the full column. See
  [the table above](@ref node_type_parameters).
- **Stagger one field**: write a one-line kernel calling [`lerp`](@ref) and dispatch
  it with [`Launcher`](@ref).
- **Differentiate**: `∂x`/`∂y` evaluated into a [`VectorField`](@ref) land
  naturally on the dual (staggered) location.
- **Compose**: `lerp`, `∂x`, `∂y`, arithmetic are ordinary Julia functions — combine
  them directly in a kernel, there is no macro to reach for.
- Boundaries are explicit (halos + `bc!`), not implicit index wrapping — see the
  note above.

=#

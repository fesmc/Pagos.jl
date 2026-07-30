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
  `ab = (Vertex, Vertex)`.
- **Stagger one field**: write a one-line kernel calling [`lerp`](@ref) and dispatch
  it with [`Launcher`](@ref).
- **Differentiate**: `∂x`/`∂y` evaluated into a [`VectorField`](@ref) land
  naturally on the dual (staggered) location.
- **Compose**: `lerp`, `∂x`, `∂y`, arithmetic are ordinary Julia functions — combine
  them directly in a kernel, there is no macro to reach for.
- Boundaries are explicit (halos + `bc!`), not implicit index wrapping — see the
  note above.

=#

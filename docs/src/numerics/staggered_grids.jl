#=

# [Staggered grids](@id staggered_grids)

Ice-dynamics solvers live on a **staggered (Arakawa C) grid**: scalars (thickness,
viscosity) sit at cell centres, while the two velocity components sit on the cell
edges, each offset by half a cell in one direction. Pagos ships a tiny,
self-contained staggered-grid toolkit in [`StaggeredGrids`](@ref) — a re-implementation
of the slice of [Oceananigans](https://github.com/CliMA/Oceananigans.jl) we actually
need: a `RectilinearGrid`, a *located* `Field`, finite-difference / interpolation
operators, and the `@at` macro to evaluate an expression at a chosen location.

The whole thing rests on two ideas:

1. **Every node carries a location** `(LX, LY, LZ)`, each entry `Center`, `Face` (or
   `Nothing` for a `Flat`, i.e. reduced, axis). A `Field` gets its location at
   construction; operators *compute* theirs from their operands.
2. **Every node is index-evaluable**: `evaluate(node, i, j, k, grid)` returns the value
   at grid point `(i, j, k)` *at that node's own location*. Interpolation is just another
   node, inserted automatically whenever an operator needs an operand at a location it
   is not yet at.

## Locations and the Pagos node naming

The familiar Pagos / Yelmo node names map directly onto a pair of horizontal
locations:

| Pagos node | location `(LX, LY)`  | lives on            | typical field            |
| ---------- | -------------------- | ------------------- | ------------------------ |
| `aa`       | `(Center, Center)`   | cell centre         | `H`, viscosity, `beta`   |
| `acx`      | `(Face,   Center)`   | east–west cell edge | `vx`                     |
| `acy`      | `(Center, Face)`     | north–south edge    | `vy`                     |
| `ab`       | `(Face,   Face)`     | cell corner         | corner-staggered `beta`  |

A horizontal-only (map-plane) grid is obtained by making the vertical axis `Flat`.

=#

using Pagos.StaggeredGrids

grid = RectilinearGrid(topology = (Periodic, Periodic, Flat),
                       size = (4, 4), x = (0, 4), y = (0, 4))

## A scalar field on aa-nodes (cell centres).
H = Field{Center, Center, Center}(grid)
H.data .= reshape(Float64.(1:16), 4, 4, 1)
H.data[:, :, 1]

#=

## Staggering an aa-field onto the ac-nodes

This is the operation you reach for constantly: a scalar known at cell centres
(`aa`) is needed on the velocity points (`acx`, `acy`). It is a half-cell average,
and you get it by asking for the field *at the target location* and materializing
with [`compute`](@ref). The lowest-level way is [`interpolate`](@ref), which inserts
the right averaging node and leaves every already-matching axis untouched:

=#

H_acx = compute(interpolate(H, (Face, Center, Center)))   # aa → acx
H_acy = compute(interpolate(H, (Center, Face, Center)))   # aa → acy

location(H_acx), location(H_acy)

#=

The result of the `acx` staggering, `H_acx[i, j] = ½ (H[i-1, j] + H[i, j])`:

=#

H_acx.data[:, :, 1]

#=

!!! note "Averaging direction (Oceananigans convention)"
    Face point `i` sits *between* centres `i-1` and `i`, so `aa → acx` averages the
    cell and its **left** neighbour: `H_acx[i,j] = ½(H[i-1,j] + H[i,j])`. The legacy
    `stagger!` used the opposite offset `½(H[i,j] + H[i+1,j])`. They differ only by the
    index convention of where face `i` lives; pick one and stay consistent. Boundaries
    are handled by the grid topology — here `Periodic`, so `i = 1` wraps to `i = 4`
    (hence `H_acx[1,1] = ½(4 + 1) = 2.5`).

## Doing it inside an expression: the `@at` macro

`interpolate` stages a single field. For a whole expression — derivatives, products,
norms — use the [`@at`](@ref) macro, which forces *every* sub-operation to evaluate at
the requested location, inserting interpolation wherever operands disagree. This is the
idiomatic way to write a staggered stencil. For example, the `x`-derivative of an
aa-field, evaluated on the `acx` velocity points:

=#

dHdx_acx = compute(@at (Face, Center, Center) ∂x(H))
dHdx_acx.data[:, :, 1]

#=

A derivative *flips* the location along its axis (`Center ↔ Face`), so `∂x` of an
aa-field naturally lands on `acx` with no extra interpolation — exactly where a
flux-divergence solver wants it. Mixed-location expressions just work; for instance an
effective strain-rate-like combination evaluated back on the centres:

=#

vx = Field{Face, Center, Center}(grid)     # an acx velocity component
vy = Field{Center, Face, Center}(grid)     # an acy velocity component
vx.data .= rand(size(vx.data)...)
vy.data .= rand(size(vy.data)...)

ε = compute(@at (Center, Center, Center) ∂x(vx)^2 + ∂y(vy)^2)
location(ε)

#=

Here `∂x(vx)` lands on `aa`, `∂y(vy)` lands on `aa`, and the sum is materialized on
`aa` — no manual staggering bookkeeping anywhere.

## In-place materialization

Both [`compute`](@ref) (allocating) and [`compute!`](@ref) (into an existing field of
the matching location) are available; use the latter in hot loops to avoid allocation:

=#

H_acx2 = Field{Face, Center, Center}(grid)
compute!(H_acx2, interpolate(H, (Face, Center, Center)))
H_acx2.data[:, :, 1]

#=

## Summary

- `aa = (Center, Center)`, `acx = (Face, Center)`, `acy = (Center, Face)`,
  `ab = (Face, Face)`; a `Flat` third axis gives a map-plane grid.
- **Stagger one field**: `compute(interpolate(field, target_location))`.
- **Stagger inside an expression**: `@at target_location expr`, then `compute`.
- Derivatives flip the location on their axis, so `∂x` of an `aa`-field lands on `acx`
  for free.
- Boundaries follow the grid `topology` (periodic wrap / bounded clamp / flat collapse);
  no halo bookkeeping is required.

See the [`StaggeredGrids`](@ref) module docstring for the full evaluation protocol and
the TODOs (stretched grids, halos, KernelAbstractions GPU port).

=#

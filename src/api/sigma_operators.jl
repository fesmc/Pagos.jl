"""
$(TYPEDSIGNATURES)

A drop-in replacement for `Chmy.∂z`, correct on a sigma-coordinate (`Chmy.FunctionAxis`)
z-axis, where `Chmy.∂z` is not. Same calling convention as Chmy's own per-dimension
derivative operators: `f` is a `Chmy.Field`, `grid` the Chmy grid it lives on, and `I` the
three grid indices at which to evaluate the derivative; the location of the result is
Chmy's own convention (opposite of `f`'s location along z, same as `f`'s along x/y).

!!! note "Why `Chmy.∂z` is wrong on a sigma axis"
    `Chmy.∂(f, grid, dim, I) = δ(f, loc, from, dim, I) * iΔ(grid, loc, dim, I)`, where
    `loc = location(f)` (the field's own location) and `from = flip(loc, dim)` (the
    location the result represents). The finite difference `δ` samples two points that
    span the *`from`*-location spacing — e.g. for `f` at `Center`, `δ = f[I] - f[I-1]`
    spans `center(I) - center(I-1)`, a `Vertex`-location distance — but Chmy scales by
    `iΔ(grid, loc, ...)`, the spacing at `loc` (`Center`), not `from`. On a `UniformAxis`,
    `spacing(ax, Center, i) == spacing(ax, Vertex, i)` identically, so the two coincide and
    the bug has no effect — which is why it went unnoticed, and why nothing needs
    correcting on the x/y axes here (always `UniformAxis`, see [`StaggeredGrid`](@ref)). On
    the sigma `FunctionAxis` the two spacings differ, and Chmy differences the right pair
    of points but divides by the wrong one — wrong by tens of percent mid-column on a
    6-layer `QuadraticSigmaTransform` grid (`test/api/sigma_operators.jl` pins the exact
    numbers). The very first/last interface (`k == 1`/`k == nz + 1`) happens to come out
    right regardless, because `_sigma_axis`'s linearly-extrapolated ghost cells make the
    two spacings coincide there by construction — which is exactly why a boundary-only
    spot check would have missed this bug entirely.

    `∂z_σ` is the one-line fix: scale by `iΔ(grid, from, dim, I)` instead of
    `iΔ(grid, loc, dim, I)`. It is built only from Chmy's exported public API
    (`location`, `flip`, `δ`, `iΔ`, `Dim`), not Chmy's unexported dispatch internals
    (`flipped`, `il`, `ir`), so a Chmy point release cannot silently break it. It reduces
    to `Chmy.∂z` bit-for-bit on any `UniformAxis` — there, `spacing(ax, Vertex, i)` and
    `spacing(ax, Center, i)` both return the same precomputed field regardless of
    location, so `loc` vs `from` truly cannot matter. On a *uniformly-spaced*
    `FunctionAxis` (e.g. `LinearSigmaTransform`) the two agree only to a few ULP, not
    bit-for-bit: `FunctionAxis` has no such fast path, so `spacing(ax, Vertex, i)` and
    `spacing(ax, Center, i)` reach the same real number by different floating-point
    routes. Either way, it is a strict correction, not a behavior change, wherever Chmy is
    already (numerically) right. Not upstreamed yet; see `roadmaps/chmy.md`, §3.

# Examples

```jldoctest
julia> layering = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 4));

julia> grid = StaggeredGrid(Float64, 4.0, 4.0, 1.0, 1.0, layering);

julia> f = Field(grid.arch, grid.grid, Center());

julia> for k in -1:(grid.nz + 2), j in -1:(grid.ny + 2), i in -1:(grid.nx + 2)
           f[i, j, k] = 3 * zcenter(grid.grid, k)   # linear in z ⟹ ∂z_σ is exactly 3 everywhere
       end

julia> round(∂z_σ(f, grid.grid, 1, 1, 1); digits = 10)   # bed interface
3.0

julia> round(∂z_σ(f, grid.grid, 1, 1, grid.nz + 1); digits = 10)   # surface interface
3.0
```
"""
@inline function ∂z_σ(f, grid, I::Vararg{Integer,3})
    loc = location(f)
    from = (loc[1], loc[2], flip(loc[3]))
    return δ(f, loc, from, Dim(3), I...) * iΔ(grid, from, Dim(3), I...)
end

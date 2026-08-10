# Shared helpers for the Chmy-native (`StaggeredGrid`/`Field`) tests.

# `Chmy.Field` allocates `2 * halo(f)` ghost cells per side (the "both halo rings" the
# docstrings below mean), and `halo(f)` is `h` (uniform on every axis) for most fields but
# a per-axis tuple for a `grid2d` field allocated `(h, h, 0)` — no ghost ring at all on `z`
# (`Pagos._halo2d`, `pagos-roadmap/memreduce.md`). A hardcoded `-1:(n + 2)` assumes the
# uniform case and writes two cells out of bounds on `z` for the tuple one; this reads the
# field's own halo per axis instead, so `k` never leaves `1:n` when `z`'s halo is `0`.
_fill_range(n, h) = (1 - 2h):(n + 2h)
_field_halo3(f) = (h = halo(f); h isa Tuple ? h : (h, h, h))

"""
    fill_analytic!(f, grid, fun)

Fill `f`'s interior *and* both halo rings from `fun(x, y)`, evaluated at `f`'s own node
class. `Chmy.set!` fills the interior only, which leaves a boundary stencil reading an
unset ghost cell; filling the halo from the analytic continuation makes the stencil
well-posed everywhere, so a test can assert over every interior point instead of only the
ones a one-sided stencil happens to reach. (In a real run those ghosts come from boundary
conditions — Phase 4 — but here the analytic continuation *is* the intended physics.)
"""
function fill_analytic!(f, grid, fun)
    loc = location(f)
    hx, hy, hz = _field_halo3(f)
    nx, ny, nz = size(interior(f))
    for k in _fill_range(nz, hz), j in _fill_range(ny, hy), i in _fill_range(nx, hx)
        x, y, _ = coord(grid, loc, i, j, k)
        f[i, j, k] = fun(x, y)
    end
    return f
end

"""
    fill_analytic3d!(f, grid, fun)

As [`fill_analytic!`](@ref) but for `fun(x, y, ζ)`, for fields that also vary with the sigma
level. Each dimension is swept over *this field's own* extent, which matters: node classes
differ in every dimension (a z-`Vertex` field has `nz + 1` layers, an x-`Vertex` field
`nx + 1` columns), so a loop range shared between components runs off the end of the
shallowest one.
"""
function fill_analytic3d!(f, grid, fun)
    loc = location(f)
    hx, hy, hz = _field_halo3(f)
    nx, ny, nz = size(interior(f))
    for k in _fill_range(nz, hz), j in _fill_range(ny, hy), i in _fill_range(nx, hx)
        x, y, ζ = coord(grid, loc, i, j, k)
        f[i, j, k] = fun(x, y, ζ)
    end
    return f
end

"""
    analytic_like(f, grid, fun)

A plain array shaped like `interior(f)`, holding `fun(x, y)` at each of `f`'s own nodes.
The reference to compare a computed field against.
"""
function analytic_like(f, grid, fun)
    out = similar(interior(f))
    loc = location(f)
    for k in axes(out, 3), j in axes(out, 2), i in axes(out, 1)
        x, y, _ = coord(grid, loc, i, j, k)
        out[i, j, k] = fun(x, y)
    end
    return out
end

"""
    convergence_rates(errors)

Observed order of accuracy between consecutive entries of `errors`, which must come from
runs on successively halved grid spacings.
"""
convergence_rates(errors) = log2.(errors[1:(end - 1)] ./ errors[2:end])

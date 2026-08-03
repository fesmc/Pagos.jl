# ---------------------------------------------------------------------------
# The generic elementwise-map kernel (`src/numerics/pointwise.jl`).
#
# `pointwise!` is shared by every pointwise physics function (`creep!`, `calving_rate!`,
# `effective_pressure!`, ...), so its overhead is paid once per such call across the whole
# package. What is worth tracking here is not the arithmetic — the caller's `f` supplies
# that — but that the launch stays allocation-free and the closure stays un-boxed: the
# kernel's own docstring records that reassigning the offset index boxes it (`Core.Box`),
# silently slower on CPU and a hard GPU compilation failure. A boxed capture would show up
# here as an allocation count, which BenchmarkTools tracks alongside the time.
#
# Reached as `Pagos.pointwise!`: it is internal machinery, not part of the exported API.
# ---------------------------------------------------------------------------

# Top-level, non-capturing: a closure defined inside `pointwise_suite` would be a different
# type on every call and drag its captures into the kernel.
_glen_viscosity(A, ε, n) = inv(2 * A^(1 / n) * ε^((n - 1) / n))

"""
    pointwise_suite(fx)

Benchmarks for `pointwise!` over 2D and 3D arrays: one array argument, and the three-array
case a real physics call site looks like.
"""
function pointwise_suite(fx)
    (; backend, T, nx, ny, nz) = fx
    g = BenchmarkGroup()

    host2d = T[BENCH_A0 * (1 + 1e-3 * (i + j)) for i in 1:nx, j in 1:ny]
    host3d = T[BENCH_A0 * (1 + 1e-3 * (i + j + k)) for i in 1:nx, j in 1:ny, k in 1:nz]

    A2 = to_backend(backend, host2d)
    e2 = to_backend(backend, fill(T(1.0e-3), nx, ny))
    o2 = to_backend(backend, zero(host2d))
    A3 = to_backend(backend, host3d)
    e3 = to_backend(backend, fill(T(1.0e-3), nx, ny, nz))
    o3 = to_backend(backend, zero(host3d))

    n = T(3)

    g["2D (2 args)"] = @benchmarkable begin
        Pagos.pointwise!(_glen_viscosity, $o2, ($A2, $e2), ($n,))
        sync!($backend)
    end
    g["3D (2 args)"] = @benchmarkable begin
        Pagos.pointwise!(_glen_viscosity, $o3, ($A3, $e3), ($n,))
        sync!($backend)
    end
    # The broadcast a hand-written call site would use instead, as the standing reference
    # for what the generic kernel costs over `@.`. Not a target to beat — `pointwise!` also
    # buys the mask/offset contract — but a gap that suddenly widens is worth seeing.
    g["3D broadcast reference"] = @benchmarkable begin
        @. $o3 = _glen_viscosity($A3, $e3, $n)
        sync!($backend)
    end

    return g
end

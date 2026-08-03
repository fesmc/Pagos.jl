# ---------------------------------------------------------------------------
# Finite-difference stencils (`src/numerics/differences.jl`).
#
# The lowest-level kernels in the package: everything above them (driving stress, velocity
# gradients, the whole PT residual) is bandwidth-bound on these, so a regression here shows
# up everywhere and is worth catching on its own.
#
# `∂x!`/`∂y!` take a *different* code path per backend (a 1-D column kernel on CPU for
# SIMD, a 2-D kernel on GPU for occupancy), so the CPU and GPU entries below are not two
# measurements of the same kernel — that is exactly why both are worth tracking.
# ---------------------------------------------------------------------------

"""
    differences_suite(fx)

Benchmarks for `∂x!`, `∂y!`, `∂x₁₂!` and both `∂x₃!` methods, on plain arrays sized from
the fixture's grid. No `Runtime` is involved: these functions read their backend off the
array itself.
"""
function differences_suite(fx)
    (; backend, T, nx, ny, nz, dx) = fx
    g = BenchmarkGroup()

    # Smooth, deterministic input — a plane wave, so no value is special (no zeros feeding
    # a division, no constant field letting a difference collapse to 0 and dodge a
    # denormal path).
    host2d = T[sinpi(2i / nx) * cospi(2j / ny) for i in 1:nx, j in 1:ny]
    host3d = T[sinpi(2i / nx) * cospi(2j / ny) * (1 + k / nz)
               for i in 1:nx, j in 1:ny, k in 1:nz]

    u2  = to_backend(backend, host2d)
    du2 = to_backend(backend, zero(host2d))
    dv2 = to_backend(backend, zero(host2d))
    u3  = to_backend(backend, host3d)
    du3 = to_backend(backend, zero(host3d))
    H   = to_backend(backend, fill(T(BENCH_H0), nx, ny))

    dxT = T(dx)
    transform = QuadraticSigmaTransform(T, nz)

    flat_x = FlatIndexing(1, nx)
    flat_y = FlatIndexing(1, ny)
    flat_z = FlatIndexing(1, nz)
    # The branchiest of the four strategies (two `mod1`s per stencil, no clamp the compiler
    # can hoist); tracked alongside the default to keep the boundary-strategy cost visible.
    periodic_x = PeriodicIndexing(1, nx)

    g["∂x!"] = @benchmarkable begin
        ∂x!($du2, $u2, $dxT, $flat_x)
        sync!($backend)
    end
    g["∂y!"] = @benchmarkable begin
        ∂y!($du2, $u2, $dxT, $flat_y)
        sync!($backend)
    end
    g["∂x! (periodic)"] = @benchmarkable begin
        ∂x!($du2, $u2, $dxT, $periodic_x)
        sync!($backend)
    end
    # Fused on GPU (one read of `u`), two sequential column kernels on CPU. The `seq/fused`
    # ratio this used to be reported as lives in `basics/differences.jl`; here it is only
    # the absolute cost of the call the physics actually makes.
    g["∂x₁₂!"] = @benchmarkable begin
        ∂x₁₂!($du2, $dv2, $u2, $dxT, $dxT, $flat_x, $flat_y)
        sync!($backend)
    end
    g["∂x₃! (uniform)"] = @benchmarkable begin
        ∂x₃!($du3, $u3, $dxT, $flat_z)
        sync!($backend)
    end
    # The sigma-coordinate method allocates a `ζ_aa` vector per call (it copies the
    # transform's midpoints onto the backend). That allocation is inside the timed region
    # on purpose: it is what a caller pays today, and if it is ever hoisted out, this
    # benchmark is where the improvement should show.
    g["∂x₃! (sigma)"] = @benchmarkable begin
        ∂x₃!($du3, $u3, $H, $transform, $flat_z)
        sync!($backend)
    end

    return g
end

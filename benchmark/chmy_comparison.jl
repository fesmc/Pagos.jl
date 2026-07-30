# Benchmark: Chmy's staggered-grid operators, with Pagos' `differences.jl` where a
# genuine equivalent still exists.
#
# `utils/oceananigans.jl` (`StaggeredGrids`) — the from-scratch Oceananigans-style
# staggered grid this benchmark used to compare against — has been deleted (see
# roadmaps/chmy.md, Phase 0): Chmy is now Pagos' staggered-grid layer, so there is
# nothing left to compare it *to* for staggered operations. What remains:
#
#   A. Directional derivatives (∂x, ∂y) — Pagos' flat, collocated `∂x₁₂!`
#      (`numerics/differences.jl`, still live and load-bearing) vs a Chmy kernel
#      calling `∂x`/`∂y` into a `VectorField`. These are genuinely different
#      operators (collocated central difference vs staggered two-point
#      difference, different output shapes) — there is no meaningful
#      cross-system agreement check, so both are validated independently
#      against the analytic derivative of a quadratic test field.
#   B. Center → Vertex-x interpolation ("stagger") via `lerp` — Chmy only,
#      validated against a linear test field (exact for linear interpolation).
#   C. A compound stencil, `∂x(H)^2 + ∂y(H)^2` evaluated back at the cell
#      centers via `lerp`, following the pattern Chmy's own test suite uses to
#      combine a `VectorField`'s components back onto their source location
#      (see Chmy's test/test_grid_operators.jl "divg" test). Chmy only,
#      validated against the analytic solution.
#
# Earlier revision bug (see roadmaps/chmy.md): case C originally read `∂x(H)`/
# `∂y(H)` directly at the *center* index without staggering back via `lerp`
# first — that silently computes a half-cell-shifted quantity. Validating
# against an analytic solution here is exactly what would have caught it.
#
# All checks compare interior points only: Chmy fields carry an explicit halo
# that nothing has filled via `bc!` in this script, so boundary-adjacent
# points read uninitialized memory by design (see docs/src/numerics/staggered_grids.jl).
#
# Run with:
#   julia -t 1    --project=. benchmark/chmy_comparison.jl   # serial CPU
#   julia -t auto --project=. benchmark/chmy_comparison.jl   # threaded CPU
#
# A CUDA comparison runs automatically iff CUDA loads and is functional.

using Pagos
using Chmy
using Chmy.Architectures
import KernelAbstractions as KA
using KernelAbstractions: @kernel, @index
using Chairmarks
using Printf

const HAS_CUDA = try
    @eval using CUDA
    CUDA.functional()
catch
    false
end

const SIZES = [(256, 256), (512, 512), (1024, 1024), (2048, 2048)]

w = 72

# ---------------------------------------------------------------------------
# Chmy kernels (mirroring Chmy's own examples/tests, e.g. examples/diffusion_2d.jl
# and test/test_grid_operators.jl).
# ---------------------------------------------------------------------------

@kernel inbounds = true function _chmy_grad!(q, H, g, O)
    I = @index(Global, NTuple)
    I = I + O
    q.x[I...] = ∂x(H, g, I...)
    q.y[I...] = ∂y(H, g, I...)
end

@kernel inbounds = true function _chmy_lerp!(out, H, g, O)
    I = @index(Global, NTuple)
    I = I + O
    out[I...] = lerp(H, location(out), g, I...)
end

@kernel inbounds = true function _chmy_grad2!(out, q, g, O)
    I = @index(Global, NTuple)
    I = I + O
    qx_c = lerp(q.x, location(out), g, I...)
    qy_c = lerp(q.y, location(out), g, I...)
    out[I...] = qx_c^2 + qy_c^2
end

# ---------------------------------------------------------------------------
# Case A: gradient (∂x and ∂y fused) — different operators, each validated
# independently against the analytic derivative of H(x, y) = x² + y².
# ---------------------------------------------------------------------------

function bench_pagos_grad(nx, ny; check = false)
    dx, dy = 1.0, 1.0
    x = (0:(nx - 1)) .* dx
    y = (0:(ny - 1)) .* dy
    u = [xi^2 + yj^2 for xi in x, yj in y]
    du1, du2 = similar(u), similar(u)
    idx1, idx2 = FlatIndexing(1, nx), FlatIndexing(1, ny)
    if check
        ∂x₁₂!(du1, du2, u, dx, dy, idx1, idx2)
        rx, ry = 2:(nx - 1), 2:(ny - 1)
        ok_x = du1[rx, ry] ≈ [2xi for xi in x[rx], _ in y[ry]]
        ok_y = du2[rx, ry] ≈ [2yj for _ in x[rx], yj in y[ry]]
        (ok_x && ok_y) || error("Pagos ∂x₁₂! disagrees with the analytic derivative")
    end
    t = (@b ∂x₁₂!($du1, $du2, $u, $dx, $dy, $idx1, $idx2)).time
    b = (@b ∂x₁₂!($du1, $du2, $u, $dx, $dy, $idx1, $idx2)).bytes
    return t, b
end

function bench_chmy_grad(backend, nx, ny; check = false)
    arch = Arch(backend)
    grid = UniformGrid(arch; origin = (0.0, 0.0), extent = (1.0, 1.0), dims = (nx, ny))
    launch = Launcher(arch, grid)
    H = Chmy.Field(arch, grid, Center())
    set!(H, grid, (x, y) -> x^2 + y^2)
    q = VectorField(arch, grid)
    if check
        launch(arch, grid, _chmy_grad! => (q, H, grid))
        xv, yv = collect(xvertices(grid)), collect(yvertices(grid))
        ic, jc = cld(nx, 2), cld(ny, 2)
        rngx, rngy = 2:(length(xv) - 1), 2:(length(yv) - 1)
        interior(q.x)[rngx, jc] ≈ 2 .* xv[rngx] ||
            error("Chmy ∂x kernel disagrees with the analytic derivative")
        interior(q.y)[ic, rngy] ≈ 2 .* yv[rngy] ||
            error("Chmy ∂y kernel disagrees with the analytic derivative")
    end
    t = (@b $launch($arch, $grid, _chmy_grad! => ($q, $H, $grid))).time
    b = (@b $launch($arch, $grid, _chmy_grad! => ($q, $H, $grid))).bytes
    return t, b
end

# ---------------------------------------------------------------------------
# Case B: Center -> Vertex-x interpolation ("stagger"), Chmy only.
# Validated against a linear field H(x, y) = x + y (exact for lerp).
# ---------------------------------------------------------------------------

function bench_chmy_lerp(backend, nx, ny; check = false)
    arch = Arch(backend)
    grid = UniformGrid(arch; origin = (0.0, 0.0), extent = (1.0, 1.0), dims = (nx, ny))
    launch = Launcher(arch, grid)
    H = Chmy.Field(arch, grid, Center())
    set!(H, grid, (x, y) -> x + y)
    out = Chmy.Field(arch, grid, (Vertex(), Center()))
    if check
        launch(arch, grid, _chmy_lerp! => (out, H, grid))
        xv, yc = collect(xvertices(grid)), collect(ycenters(grid))
        jc = cld(ny, 2)
        rng = 2:(length(xv) - 1)
        interior(out)[rng, jc] ≈ xv[rng] .+ yc[jc] ||
            error("Chmy lerp kernel disagrees with the analytic (linear) solution")
    end
    t = (@b $launch($arch, $grid, _chmy_lerp! => ($out, $H, $grid))).time
    b = (@b $launch($arch, $grid, _chmy_lerp! => ($out, $H, $grid))).bytes
    return t, b
end

# ---------------------------------------------------------------------------
# Case C: compound stencil ∂x(H)^2 + ∂y(H)^2, staggered back onto the centers
# it started from — the location-faithful pattern from
# docs/src/numerics/staggered_grids.jl. Chmy only. For H(x, y) = x² + y², the
# analytic result has the closed form 4·H, so the check reuses H directly.
# ---------------------------------------------------------------------------

function bench_chmy_grad2(backend, nx, ny; check = false)
    arch = Arch(backend)
    grid = UniformGrid(arch; origin = (0.0, 0.0), extent = (1.0, 1.0), dims = (nx, ny))
    launch = Launcher(arch, grid)
    H = Chmy.Field(arch, grid, Center())
    set!(H, grid, (x, y) -> x^2 + y^2)
    q = VectorField(arch, grid)
    launch(arch, grid, _chmy_grad! => (q, H, grid))
    out = Chmy.Field(arch, grid, Center())
    if check
        launch(arch, grid, _chmy_grad2! => (out, q, grid))
        r = 2:(nx - 1)
        jc = cld(ny, 2)
        interior(out)[r, jc] ≈ 4 .* interior(H)[r, jc] ||
            error("Chmy compound-stencil kernel disagrees with the analytic solution")
    end
    t = (@b $launch($arch, $grid, _chmy_grad2! => ($out, $q, $grid))).time
    b = (@b $launch($arch, $grid, _chmy_grad2! => ($out, $q, $grid))).bytes
    return t, b
end

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

function print_title(title)
    println("=" ^ w)
    println(" " * title)
    println("=" ^ w)
end

function run_case(name, bench_fn, backend)
    print_title(name)
    @printf("  %-14s %14s %14s\n", "Grid", "time [us]", "bytes")
    println("-" ^ w)
    for (nx, ny) in SIZES
        t, b = bench_fn(backend, nx, ny; check = true)
        @printf("  %-14s %14.2f %14d\n", "$(nx)×$(ny)", t * 1e6, b)
    end
    println("=" ^ w)
    println()
end

function run_pair(name, pagos_fn, chmy_fn, backend)
    print_title(name)
    @printf("  %-14s %14s %14s %14s %14s %9s\n",
        "Grid", "Pagos [us]", "Pagos [B]", "Chmy [us]", "Chmy [B]", "ratio")
    println("-" ^ w)
    for (nx, ny) in SIZES
        tp, bp = pagos_fn(nx, ny; check = true)
        tc, bc = chmy_fn(backend, nx, ny; check = true)
        @printf("  %-14s %14.2f %14d %14.2f %14d %8.2fx\n",
            "$(nx)×$(ny)", tp * 1e6, bp, tc * 1e6, bc, tp / tc)
    end
    println("=" ^ w)
    println()
end

println("Julia threads = ", Threads.nthreads())
println("KernelAbstractions CPU backend\n")

run_pair("A. gradient: Pagos ∂x₁₂! vs Chmy ∂x/∂y (different operators, both validated)",
    bench_pagos_grad, bench_chmy_grad, KA.CPU())
run_case("B. stagger: Chmy lerp", bench_chmy_lerp, KA.CPU())
run_case("C. compound: ∂x(H)^2+∂y(H)^2 staggered back to centers, Chmy", bench_chmy_grad2, KA.CPU())

if HAS_CUDA
    println("CUDA device = ", CUDA.name(CUDA.device()), "\n")
    run_pair("A. gradient [GPU]", bench_pagos_grad, bench_chmy_grad, CUDABackend())
    run_case("B. stagger [GPU]", bench_chmy_lerp, CUDABackend())
    run_case("C. compound [GPU]", bench_chmy_grad2, CUDABackend())
else
    println("CUDA not available / not functional — GPU comparison skipped.")
end

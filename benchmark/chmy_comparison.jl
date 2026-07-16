# Benchmark comparison: Chmy's staggered-grid operators vs Pagos' own
# differential-operator implementations.
#
# Pagos currently has two, unrelated operator layers:
#   1. `numerics/differences.jl` — flat, collocated central differences on a
#      plain array (∂x!, ∂y!, ∂x₁₂!), boundary behaviour controlled by
#      `AbstractIndexing`. Output stays the same size/location as the input.
#   2. `utils/oceananigans.jl` (`StaggeredGrids`) — a small from-scratch
#      Oceananigans-style staggered grid: `Field`s carry a `(Center/Face)`
#      location, `∂x`/`∂y`/interpolation are lazy nodes, and `@at` inserts the
#      interpolation needed to reconcile locations automatically.
#
# Chmy (already a declared dependency, see Project.toml) provides a
# `KernelAbstractions`-based, multi-GPU/MPI-capable staggered grid with the
# same conceptual pieces (`Center()`/`Vertex()` locations, `∂x`/`∂y`, `lerp`).
# Unlike `StaggeredGrids`, Chmy's operators are *not* lazy/auto-interpolating:
# `∂x(f, grid, I...)` just evaluates the stencil at whatever index `I` you
# hand it — it is on you to `lerp` mismatched locations back together.
#
# Three matched pairs are benchmarked here:
#   A. Directional derivatives (∂x, ∂y) fused into a gradient — Pagos' flat
#      `∂x₁₂!` vs a Chmy kernel calling `∂x`/`∂y` into a `VectorField`.
#   B. Center → Vertex-x interpolation ("stagger") — Pagos `StaggeredGrids`
#      `interpolate` vs Chmy `lerp`.
#   C. A compound stencil, `∂x(H)^2 + ∂y(H)^2` evaluated back at the
#      cell centers — Pagos `StaggeredGrids` `@at` vs a hand-written Chmy
#      kernel (no `@at` equivalent exists in Chmy: the user writes the
#      combination directly).
#
# Run with:
#   julia -t 1    --project=. benchmark/chmy_comparison.jl   # serial CPU
#   julia -t auto --project=. benchmark/chmy_comparison.jl   # threaded CPU
#
# A CUDA comparison runs automatically iff CUDA loads and is functional.

using Pagos
import Pagos.StaggeredGrids as SG
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

w = 96

# ---------------------------------------------------------------------------
# Chmy kernels (no ready-made "gradient"/"stagger" API call — these one-liners
# are the idiomatic way to use Chmy's operators, mirroring Chmy's own examples
# such as examples/diffusion_2d.jl and test/test_grid_operators.jl).
# ---------------------------------------------------------------------------

@kernel inbounds = true function _chmy_grad!(qx, qy, C, g, O)
    I = @index(Global, NTuple)
    I = I + O
    qx[I...] = ∂x(C, g, I...)
    qy[I...] = ∂y(C, g, I...)
end

@kernel inbounds = true function _chmy_lerp!(out, H, g, O)
    I = @index(Global, NTuple)
    I = I + O
    out[I...] = lerp(H, (Vertex(), Center()), g, I...)
end

@kernel inbounds = true function _chmy_grad2!(out, H, g, O)
    I = @index(Global, NTuple)
    I = I + O
    out[I...] = ∂x(H, g, I...)^2 + ∂y(H, g, I...)^2
end

# ---------------------------------------------------------------------------
# Case A: gradient (∂x and ∂y fused)
# ---------------------------------------------------------------------------

function bench_pagos_grad(nx, ny)
    u = rand(nx, ny)
    du1, du2 = similar(u), similar(u)
    idx1, idx2 = FlatIndexing(1, nx), FlatIndexing(1, ny)
    t = (@b ∂x₁₂!($du1, $du2, $u, 1.0, 1.0, $idx1, $idx2)).time
    b = (@b ∂x₁₂!($du1, $du2, $u, 1.0, 1.0, $idx1, $idx2)).bytes
    return t, b
end

function bench_chmy_grad(backend, nx, ny)
    arch = Arch(backend)
    grid = UniformGrid(arch; origin = (0.0, 0.0), extent = (1.0, 1.0), dims = (nx, ny))
    launch = Launcher(arch, grid)
    C = Chmy.Field(arch, grid, Center())
    set!(C, grid, (x, y) -> x + y)
    qx = Chmy.Field(arch, grid, (Vertex(), Center()))
    qy = Chmy.Field(arch, grid, (Center(), Vertex()))
    t = (@b $launch($arch, $grid, _chmy_grad! => ($qx, $qy, $C, $grid))).time
    b = (@b $launch($arch, $grid, _chmy_grad! => ($qx, $qy, $C, $grid))).bytes
    return t, b
end

# ---------------------------------------------------------------------------
# Case B: Center -> Vertex-x interpolation ("stagger")
# ---------------------------------------------------------------------------

function bench_pagos_lerp(nx, ny)
    grid = SG.RectilinearGrid(topology = (SG.Periodic, SG.Periodic, SG.Flat),
                               size = (nx, ny), x = (0, 1), y = (0, 1))
    H = SG.Field{SG.Center, SG.Center, SG.Center}(grid)
    H.data .= rand(nx, ny, 1)
    op = SG.interpolate(H, (SG.Face, SG.Center, SG.Center))
    out = SG.compute(op)
    t = (@b SG.compute!($out, $op)).time
    b = (@b SG.compute!($out, $op)).bytes
    return t, b
end

function bench_chmy_lerp(backend, nx, ny)
    arch = Arch(backend)
    grid = UniformGrid(arch; origin = (0.0, 0.0), extent = (1.0, 1.0), dims = (nx, ny))
    launch = Launcher(arch, grid)
    H = Chmy.Field(arch, grid, Center())
    set!(H, grid, (x, y) -> x + y)
    out = Chmy.Field(arch, grid, (Vertex(), Center()))
    t = (@b $launch($arch, $grid, _chmy_lerp! => ($out, $H, $grid))).time
    b = (@b $launch($arch, $grid, _chmy_lerp! => ($out, $H, $grid))).bytes
    return t, b
end

# ---------------------------------------------------------------------------
# Case C: compound stencil ∂x(H)^2 + ∂y(H)^2
# ---------------------------------------------------------------------------

function bench_pagos_grad2(nx, ny)
    grid = SG.RectilinearGrid(topology = (SG.Periodic, SG.Periodic, SG.Flat),
                               size = (nx, ny), x = (0, 1), y = (0, 1))
    H = SG.Field{SG.Center, SG.Center, SG.Center}(grid)
    H.data .= rand(nx, ny, 1)
    op = SG.@at (SG.Center, SG.Center, SG.Center) ∂x(H)^2 + ∂y(H)^2
    out = SG.compute(op)
    t = (@b SG.compute!($out, $op)).time
    b = (@b SG.compute!($out, $op)).bytes
    return t, b
end

function bench_chmy_grad2(backend, nx, ny)
    arch = Arch(backend)
    grid = UniformGrid(arch; origin = (0.0, 0.0), extent = (1.0, 1.0), dims = (nx, ny))
    launch = Launcher(arch, grid)
    H = Chmy.Field(arch, grid, Center())
    set!(H, grid, (x, y) -> x + y)
    out = Chmy.Field(arch, grid, Center())
    t = (@b $launch($arch, $grid, _chmy_grad2! => ($out, $H, $grid))).time
    b = (@b $launch($arch, $grid, _chmy_grad2! => ($out, $H, $grid))).bytes
    return t, b
end

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

function print_header(title)
    println("=" ^ w)
    println(" " * title)
    println("=" ^ w)
    @printf("  %-14s %14s %14s %14s %14s %10s\n",
        "Grid", "Pagos [us]", "Pagos [B]", "Chmy [us]", "Chmy [B]", "speedup")
    println("-" ^ w)
end

function run_case(name, pagos_fn, chmy_fn, backend)
    print_header(name)
    for (nx, ny) in SIZES
        tp, bp = pagos_fn(nx, ny)
        tc, bc = chmy_fn(backend, nx, ny)
        @printf("  %-14s %14.2f %14d %14.2f %14d %9.2fx\n",
            "$(nx)×$(ny)", tp * 1e6, bp, tc * 1e6, bc, tp / tc)
    end
    println("=" ^ w)
    println()
end

println("Julia threads = ", Threads.nthreads())
println("KernelAbstractions CPU backend\n")

run_case("A. gradient: Pagos ∂x₁₂! vs Chmy ∂x/∂y", bench_pagos_grad, bench_chmy_grad, KA.CPU())
run_case("B. stagger: Pagos StaggeredGrids.interpolate vs Chmy lerp", bench_pagos_lerp, bench_chmy_lerp, KA.CPU())
run_case("C. compound: Pagos @at ∂x(H)^2+∂y(H)^2 vs Chmy hand-written kernel", bench_pagos_grad2, bench_chmy_grad2, KA.CPU())

if HAS_CUDA
    println("CUDA device = ", CUDA.name(CUDA.device()), "\n")
    run_case("A. gradient [GPU]", bench_pagos_grad, bench_chmy_grad, CUDABackend())
    run_case("B. stagger [GPU]", bench_pagos_lerp, bench_chmy_lerp, CUDABackend())
    run_case("C. compound [GPU]", bench_pagos_grad2, bench_chmy_grad2, CUDABackend())
else
    println("CUDA not available / not functional — GPU comparison skipped.")
end

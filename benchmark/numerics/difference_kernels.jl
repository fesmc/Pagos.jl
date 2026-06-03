using Pagos
using Chairmarks
using Printf
using CUDA

# ---------------------------------------------------------------------------
# Benchmark: ∂x₁ / ∂x₂ / ∂x₁₂ across CPU (serial & threaded) and CUDA GPU
#
# Thread count is fixed at Julia startup.  Run with:
#   julia -t 1  benchmark/numerics/difference_kernels.jl   # serial CPU
#   julia -t 4  benchmark/numerics/difference_kernels.jl   # 4-thread CPU
#
# The functions call KernelAbstractions.synchronize internally, so the
# benchmark captures the full round-trip time including synchronization.
# ---------------------------------------------------------------------------

const HAS_CUDA = CUDA.functional()
const SIZES    = [(128, 128), (256, 256), (512, 512), (1024, 1024)]

w = 80

function bench_cpu(nx, ny)
    u   = rand(Float64, nx, ny)
    du₁ = similar(u)
    du₂ = similar(u)
    idx₁ = FlatIndexing(1, nx)
    idx₂ = FlatIndexing(1, ny)
    t1  = (@b ∂x₁!($du₁, $u, 1.0, $idx₁)).time
    t2  = (@b ∂x₂!($du₂, $u, 1.0, $idx₂)).time
    t12 = (@b ∂x₁₂!($du₁, $du₂, $u, 1.0, 1.0, $idx₁, $idx₂)).time
    return t1, t2, t12
end

function bench_cuda(nx, ny)
    u   = CUDA.rand(Float64, nx, ny)
    du₁ = similar(u)
    du₂ = similar(u)
    idx₁ = FlatIndexing(1, nx)
    idx₂ = FlatIndexing(1, ny)
    t1  = (@b ∂x₁!($du₁, $u, 1.0, $idx₁)).time
    t2  = (@b ∂x₂!($du₂, $u, 1.0, $idx₂)).time
    t12 = (@b ∂x₁₂!($du₁, $du₂, $u, 1.0, 1.0, $idx₁, $idx₂)).time
    return t1, t2, t12
end

function print_header(title)
    println("=" ^ w)
    println(" " * title)
    println("=" ^ w)
    @printf("  %-14s %11s %11s %11s %9s\n",
        "Grid", "∂x₁ [μs]", "∂x₂ [μs]", "∂x₁₂ [μs]", "seq/fused")
    println("-" ^ w)
end

function print_row(nx, ny, t1, t2, t12)
    @printf("  %-14s %11.1f %11.1f %11.1f %9.2fx\n",
        "$(nx)×$(ny)",
        t1  * 1e6,
        t2  * 1e6,
        t12 * 1e6,
        (t1 + t2) / t12)
end

# ---------------------------------------------------------------------------
# CPU benchmark
# ---------------------------------------------------------------------------

ncpu = Threads.nthreads()
print_header("∂x₁ / ∂x₂ / ∂x₁₂ — CPU ($ncpu thread$(ncpu == 1 ? "" : "s"), KernelAbstractions CPU backend)")

for (nx, ny) in SIZES
    t1, t2, t12 = bench_cpu(nx, ny)
    print_row(nx, ny, t1, t2, t12)
end

println("=" ^ w)
println()

# ---------------------------------------------------------------------------
# CUDA GPU benchmark
# ---------------------------------------------------------------------------

if HAS_CUDA
    devname = CUDA.name(CUDA.device())
    print_header("∂x₁ / ∂x₂ / ∂x₁₂ — CUDA GPU ($devname)")

    for (nx, ny) in SIZES
        t1, t2, t12 = bench_cuda(nx, ny)
        print_row(nx, ny, t1, t2, t12)
    end

    println("=" ^ w)
    println()
else
    println("GPU benchmark skipped: no CUDA-capable device found (CUDA.functional() == false)")
    println()
end

println("Notes:")
println("  Run with `julia -t 1` for serial and `julia -t 4` for 4-threaded CPU results.")
println("  ∂x₁₂ reads u once; ∂x₁ + ∂x₂ reads u twice — fused wins when u > L3 cache.")
println("  Synchronization is included in all timings (KernelAbstractions.synchronize).")

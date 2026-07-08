# Benchmark for the StaggeredGrids utils (src/utils/oceananigans.jl).
#
# Two questions:
#   1. Is the serial `compute!` allocation-free in its hot loop? (it should report 0 B)
#   2. How does a parallel KernelAbstractions kernel over the same `evaluate` tree
#      compare to the serial CPU loop?
#
# `compute!` in the module is a plain serial `@inbounds` loop. The parallel variant
# below (`compute_ka!`) wraps the *same* pure, type-stable `evaluate` tree in a
# KernelAbstractions kernel — exactly the GPU-ready port noted in the module TODOs.
# On the CPU() backend it multithreads across `Threads.nthreads()`, so run with
# several threads to see scaling:
#
#   julia -t auto --project=. benchmark/oceananigans.jl
#
# A GPU comparison runs automatically *iff* CUDA loads and `CUDA.functional()`. CUDA
# is only a weakdep of Pagos, so run the GPU case from an environment that has both
# Pagos and CUDA (e.g. the benchmark project with Pagos dev-added); otherwise it is
# silently skipped and only the CPU serial-vs-parallel comparison runs.

using Pagos
using Pagos.StaggeredGrids
import KernelAbstractions as KA
using KernelAbstractions: @kernel, @index
using Adapt
using Printf

# Load CUDA only if present; gate the whole GPU section on a functional device.
const HAS_CUDA = try
    @eval using CUDA
    CUDA.functional()
catch
    false
end

# --- Parallel compute!: one KA kernel over the (pure, type-stable) evaluate tree ---
@kernel function _compute_kernel!(out, op, grid)
    i, j, k = @index(Global, NTuple)
    @inbounds out[i, j, k] = evaluate(op, i, j, k, grid)
end

function compute_ka!(out::Field, op, backend = KA.CPU())
    grid = out.grid
    _compute_kernel!(backend)(out.data, op, grid; ndrange = size(out.data))
    KA.synchronize(backend)
    return out
end

# min-of-N timer (same idea as benchmark/stress.jl; kernels are ms-scale)
function bench(f, args...; samples = 200)
    f(args...)                       # warmup / compile
    best = Inf
    for _ in 1:samples
        best = min(best, @elapsed f(args...))
    end
    return best
end

# allocations of a single warm call
function allocs(f, args...)
    f(args...)                       # warmup / compile
    return @allocated f(args...)
end

# Measure one concrete operation (kept in its own function so `op` is type-stable here).
function measure(label, op, grid, backend)
    out_s = compute(op)                        # serial result + warmup, allocates the buffer
    out_p = Field{location(op)...}(grid)       # matching-location buffer for the parallel run
    compute_ka!(out_p, op, backend)            # warmup
    agree = out_s.data ≈ out_p.data            # the two paths must give the same field

    a_s = allocs(compute!,    out_s, op)
    a_p = allocs(compute_ka!, out_p, op, backend)
    t_s = bench(compute!,    out_s, op)
    t_p = bench(compute_ka!, out_p, op, backend)

    @printf("  %-18s N=%9d   serial %8.3f ms /%7d B    parallel %8.3f ms /%8d B    speedup %5.2fx   %s\n",
            label, length(out_s.data), 1e3t_s, a_s, 1e3t_p, a_p, t_s / t_p,
            agree ? "match" : "MISMATCH!")
    return nothing
end

function run_case(n; backend = KA.CPU())
    grid = RectilinearGrid(topology = (Periodic, Periodic, Flat),
                           size = (n, n), x = (0, 1), y = (0, 1))
    H = Field{Center, Center, Center}(grid)
    H.data .= rand(eltype(grid), size(H.data)...)

    println("grid $(n)x$(n):")
    # 1. The staggering the docs showcase: aa-node scalar onto the acx velocity points.
    measure("aa->acx stagger", interpolate(H, (Face, Center, Center)), grid, backend)
    # 2. A heavier compound stencil (two derivatives + interpolation back to centres).
    measure("dx^2+dy^2 @ aa", (@at (Center, Center, Center) ∂x(H)^2 + ∂y(H)^2), grid, backend)
    return nothing
end

# Optional GPU case: serial CPU baseline vs the same kernel on a CUDA device.
# The op tree is moved to the device with `adapt(CuArray, op)` (see the Adapt rules
# in src/utils/oceananigans.jl), so the leaf `Field` data becomes a `CuArray`.
function run_case_gpu(n)
    backend = CUDABackend()
    grid = RectilinearGrid(topology = (Periodic, Periodic, Flat),
                           size = (n, n), x = (0, 1), y = (0, 1))
    H = Field{Center, Center, Center}(grid)
    H.data .= rand(eltype(grid), size(H.data)...)

    println("grid $(n)x$(n) [GPU]:")
    measure_gpu("aa->acx stagger", interpolate(H, (Face, Center, Center)), grid, backend)
    measure_gpu("dx^2+dy^2 @ aa", (@at (Center, Center, Center) ∂x(H)^2 + ∂y(H)^2), grid, backend)
    return nothing
end

function measure_gpu(label, op_cpu, grid, backend)
    out_ref = compute(op_cpu)                                     # CPU reference + serial buffer
    op  = adapt(CuArray, op_cpu)                                  # move leaf data onto the device
    out = Field{location(op_cpu)...}(grid, CuArray(zero(out_ref.data)))
    compute_ka!(out, op, backend)                                 # warmup / kernel compile
    agree = Array(out.data) ≈ out_ref.data

    t_cpu = bench(compute!,    out_ref, op_cpu)                   # serial CPU baseline
    t_gpu = bench(compute_ka!, out, op, backend)                 # GPU (compute_ka! synchronizes)
    a_gpu = allocs(compute_ka!, out, op, backend)                # host-side launch allocations

    @printf("  %-18s N=%9d   serial CPU %8.3f ms    GPU %8.3f ms /%7d B(host)    speedup %7.1fx   %s\n",
            label, length(out_ref.data), 1e3t_cpu, 1e3t_gpu, a_gpu, t_cpu / t_gpu,
            agree ? "match" : "MISMATCH!")
    return nothing
end

println("Julia threads = ", Threads.nthreads())
println("parallel backend = KernelAbstractions ", KA.CPU(), "\n")

run_case(512)
run_case(1024)
run_case(2048)

if HAS_CUDA
    println("\nCUDA device = ", CUDA.name(CUDA.device()), "\n")
    run_case_gpu(512)
    run_case_gpu(1024)
    run_case_gpu(2048)
    run_case_gpu(4096)
else
    println("\nCUDA not available / not functional — GPU benchmark skipped.")
end

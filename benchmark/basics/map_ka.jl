# ---------------------------------------------------------------------------
# Benchmark: map! vs KernelAbstractions
# Task: regularized Coulomb basal friction β over a 512×512 grid.
#
#   β = c_bed · (‖v‖ / (‖v‖ + v₀))^q / ‖v‖,   ‖v‖ = √(vx²+vy²) + v_reg
#
# Run:
#   julia -t 1  benchmark/numerics/map_ka.jl   # serial CPU
#   julia -t 4  benchmark/numerics/map_ka.jl   # 4-thread CPU
#
# map!    → single-threaded on CPU; dispatches to a CuArray kernel on GPU
# KA      → Threads.nthreads() on CPU; fused kernel on GPU
#
# GPU timings use CUDA.@elapsed (GPU-side events), which correctly captures
# async kernel execution time.  CPU timings use Chairmarks @b (wall clock).
# ---------------------------------------------------------------------------

using Pkg
Pkg.activate(".")
using Chairmarks
using CUDA
using Printf
using KernelAbstractions

const HAS_CUDA = try
    using CUDA
    CUDA.functional()
catch
    false
end

const N    = 512
const T    = Float64
const V0   = T(1e2)    # transition velocity [m/yr]
const VREG = T(1e-3)   # regularization velocity [m/yr]
const Q    = T(0.2)    # Coulomb sliding exponent

# ── scalar formula ──────────────────────────────────────────────────────────
@inline function coulomb_beta(c, vx, vy, v0, vreg, q)
    vn = sqrt(vx * vx + vy * vy) + vreg
    return c * (vn / (vn + v0))^q / vn
end

# ── KernelAbstractions kernel ───────────────────────────────────────────────
@kernel function _coulomb_ka!(β, @Const(c_bed), @Const(vx), @Const(vy), v0, vreg, q)
    i, j = @index(Global, NTuple)
    β[i, j] = coulomb_beta(c_bed[i, j], vx[i, j], vy[i, j], v0, vreg, q)
end

# ── launchers ───────────────────────────────────────────────────────────────
function run_map!(β, c, vx, vy)
    map!((ci, vxi, vyi) -> coulomb_beta(ci, vxi, vyi, V0, VREG, Q), β, c, vx, vy)
end

function run_ka!(β, c, vx, vy, backend)
    kernel = _coulomb_ka!(backend, (16, 16))
    kernel(β, c, vx, vy, V0, VREG, Q; ndrange = size(β))
    KernelAbstractions.synchronize(backend)
end

# ── GPU timing: CUDA.@elapsed measures GPU-side event time, not wall clock ──
function gpu_min_time(f; warmup = 5, N = 200)
    for _ in 1:warmup
        f()
    end
    CUDA.synchronize()
    t = Inf
    for _ in 1:N
        t = min(t, CUDA.@elapsed f())
    end
    return t
end

# ── reporting ───────────────────────────────────────────────────────────────
const W = 66

function print_header(title)
    println("=" ^ W)
    println(" " * title)
    println("=" ^ W)
    @printf("  %-18s %12s %10s\n", "Method", "time [μs]", "vs map!")
    println("-" ^ W)
end

function print_row(name, t_method, t_ref)
    @printf("  %-18s %12.1f %10s\n", name, t_method * 1e6,
        @sprintf("%.2f×", t_ref / t_method))
end

# ── CPU ─────────────────────────────────────────────────────────────────────
let
    c  = rand(T, N, N)
    vx = rand(T, N, N) .* T(500)
    vy = rand(T, N, N) .* T(500)
    β  = similar(c)
    backend = KernelAbstractions.CPU()
    ncpu = Threads.nthreads()

    print_header("Coulomb β — CPU ($ncpu thread$(ncpu == 1 ? "" : "s"), $(N)×$(N) $T)")

    run_map!(β, c, vx, vy); run_ka!(β, c, vx, vy, backend)

    t_map = (@b run_map!($β, $c, $vx, $vy)).time
    t_ka  = (@b run_ka!($β, $c, $vx, $vy, $backend)).time

    print_row("map!", t_map, t_map)
    print_row("KA",   t_ka,  t_map)
    println("=" ^ W); println()
end

# ── GPU ─────────────────────────────────────────────────────────────────────
if HAS_CUDA
    let
        devname = CUDA.name(CUDA.device())
        c  = CUDA.rand(T, N, N)
        vx = CUDA.rand(T, N, N) .* T(500)
        vy = CUDA.rand(T, N, N) .* T(500)
        β  = similar(c)
        backend = CUDABackend()

        print_header("Coulomb β — GPU ($devname, $(N)×$(N) $T)")

        t_map = gpu_min_time(() -> run_map!(β, c, vx, vy))
        t_ka  = gpu_min_time(() -> run_ka!(β, c, vx, vy, backend))

        print_row("map!", t_map, t_map)
        print_row("KA",   t_ka,  t_map)
        println("=" ^ W); println()
    end
else
    println("GPU benchmark skipped: no CUDA-capable device (CUDA.functional() == false).")
end

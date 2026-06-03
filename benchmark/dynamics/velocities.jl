using Pagos
using BenchmarkTools
using LinearSolve, SparseArrays
using Printf
using CUDA, CUDA.CUSPARSE

const HAS_CUDA = CUDA.functional()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function cpu_inputs(nx, ny; T = Float64)
    return (
        rand(T, nx, ny),   # N
        rand(T, nx, ny),   # N_ab
        rand(T, nx, ny),   # ux
        rand(T, nx, ny),   # uy
        rand(T, nx, ny),   # taud_acx
        rand(T, nx, ny),   # taud_acy
        rand(T, nx, ny),   # β_acx
        rand(T, nx, ny),   # β_acy
    )
end

function gpu_inputs(cpu_ins)
    return CuArray.(cpu_ins)
end

function gpu_solver(lsd2::LinearSolver2D_v2)
    return LinearSolver2D_v2(
        lsd2.dynamics,
        lsd2.resolution_params,
        CuArray(lsd2.u),
        CuArray(lsd2.u0),
        CuArray(lsd2.b),
        CuSparseMatrixCSC(lsd2.A),
        CuArray(lsd2.perm),
        lsd2.i_idx,
        lsd2.j_idx,
    )
end

# ---------------------------------------------------------------------------
# Benchmark loop
# ---------------------------------------------------------------------------
# For each grid size we measure:
#   v1  : populate_vectors! (COO fill) + sparse(Ai,Aj,Av) assembly
#   v2  : populate_vectors! on CPU via KA CPU() kernel, no assembly step
#   gpu : populate_vectors! on GPU via KA CUDABackend() kernel
# ---------------------------------------------------------------------------

const SIZES = [(50, 50), (128, 128), (256, 256), (512, 512)]
const T     = Float64

w = 72
println("=" ^ w)
println(" populate_vectors! benchmark — v1 (COO) vs v2 CPU (CSC) vs v2 GPU")
println("=" ^ w)
@printf "  %-18s %11s %11s %11s %9s %9s\n" \
    "Grid (2×n DOF)" "v1 [μs]" "v2-CPU [μs]" "v2-GPU [μs]" "v2/v1" "GPU/CPU"
println("-" ^ w)

for (nx, ny) in SIZES
    rp = ResolutionParameters(2, nx, ny, 1.0, 1.0)
    ins = cpu_inputs(nx, ny; T)
    N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy = ins

    # ---- v1 -----------------------------------------------------------------
    lsd1 = LinearSolver2D(rp; T)
    populate_vectors!(lsd1, N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy, DIVA())
    t_fill  = @belapsed populate_vectors!(
        $lsd1, $N, $N_ab, $ux, $uy, $taud_acx, $taud_acy, $β_acx, $β_acy, DIVA()
    )
    t_assem = @belapsed sparse($lsd1.Ai, $lsd1.Aj, $lsd1.Av)
    t_v1    = t_fill + t_assem

    # ---- v2 CPU -------------------------------------------------------------
    lsd2 = LinearSolver2D_v2(rp, DIVA(); T)
    t_v2_cpu = @belapsed populate_vectors!(
        $lsd2, $N, $N_ab, $ux, $uy, $taud_acx, $taud_acy, $β_acx, $β_acy
    )

    # ---- v2 GPU -------------------------------------------------------------
    if HAS_CUDA
        ins_cu   = gpu_inputs(ins)
        N_cu, N_ab_cu, ux_cu, uy_cu, taud_acx_cu, taud_acy_cu, β_acx_cu, β_acy_cu = ins_cu
        lsd2_gpu = gpu_solver(lsd2)

        # one warmup pass to trigger JIT compilation
        populate_vectors!(lsd2_gpu, N_cu, N_ab_cu, ux_cu, uy_cu,
                          taud_acx_cu, taud_acy_cu, β_acx_cu, β_acy_cu)
        CUDA.synchronize()

        t_v2_gpu = @belapsed begin
            populate_vectors!(
                $lsd2_gpu, $N_cu, $N_ab_cu, $ux_cu, $uy_cu,
                $taud_acx_cu, $taud_acy_cu, $β_acx_cu, $β_acy_cu,
            )
            CUDA.synchronize()
        end

        @printf "  %-18s %11.1f %11.1f %11.1f %8.1fx %8.1fx\n" \
            "$(nx)×$(ny) ($(2*nx*ny))" \
            t_v1 * 1e6 \
            t_v2_cpu * 1e6 \
            t_v2_gpu * 1e6 \
            t_v1 / t_v2_cpu \
            t_v2_cpu / t_v2_gpu
    else
        @printf "  %-18s %11.1f %11.1f %11s %8.1fx %9s\n" \
            "$(nx)×$(ny) ($(2*nx*ny))" \
            t_v1 * 1e6 \
            t_v2_cpu * 1e6 \
            "N/A" \
            t_v1 / t_v2_cpu \
            "N/A"
    end
end

println("=" ^ w)
println()
println("Notes:")
println("  v1      COO fill → sparse(Ai,Aj,Av) reassembly on every solve")
println("  v2-CPU  KA CPU() kernel → direct nzval write, no reassembly")
println("  v2-GPU  KA CUDABackend() kernel, CuSparseMatrixCSC; includes CUDA.synchronize()")
HAS_CUDA || println("  GPU benchmarks skipped: no CUDA-capable device found (CUDA.functional() == false)")

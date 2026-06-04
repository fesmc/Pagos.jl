using Pagos
using Chairmarks
using SparseArrays
using Printf
using CUDA, CUDA.CUSPARSE
using CUDSS

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

function gpu_solver(lsd2::LinearDynamicsSolver2D)
    return LinearDynamicsSolver2D(
        lsd2.dynamics,
        lsd2.resolution_params,
        CuArray(lsd2.u),
        CuArray(lsd2.u0),
        CuArray(lsd2.b),
        CuSparseMatrixCSC(lsd2.A),
        CuArray(lsd2.perm),
        lsd2.i_idx,
        lsd2.j_idx,
        Ref{Any}(nothing),
    )
end

# CUDSS requires CSR format. The `perm` array maps COO fill-order → CSC nzval positions;
# we remap it to CSR positions by sorting the CSC entries by (row, col).
function gpu_solver_csr(lsd2::LinearDynamicsSolver2D)
    A_cpu = lsd2.A                              # SparseMatrixCSC on CPU
    Ai_csc, Aj_csc, _ = findnz(A_cpu)          # (row, col) in CSC nzval order
    csc_to_csr = invperm(sortperm(collect(zip(Ai_csc, Aj_csc))))
    perm_csr   = csc_to_csr[Array(lsd2.perm)]  # remap perm to CSR nzval positions
    return LinearDynamicsSolver2D(
        lsd2.dynamics,
        lsd2.resolution_params,
        CuArray(lsd2.u),
        CuArray(lsd2.u0),
        CuArray(lsd2.b),
        CuSparseMatrixCSR(A_cpu),
        CuArray(perm_csr),
        lsd2.i_idx,
        lsd2.j_idx,
        Ref{Any}(nothing),
    )
end

# ---------------------------------------------------------------------------
# Setup: build all solver instances and inputs upfront
# ---------------------------------------------------------------------------

const SIZES = [(50, 50), (128, 128), (256, 256), (512, 512)]
const T = Float64

all_data = map(SIZES) do (nx, ny)
    rp  = ResolutionParameters(2, nx, ny, 1.0, 1.0)
    ins = cpu_inputs(nx, ny; T)
    N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy = ins

    lsd1 = LegacyLinearDynamicsSolver2D(rp; T)
    populate_vectors!(lsd1, N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy, DIVADynamicsXY())

    lsd2 = LinearDynamicsSolver2D(rp, DIVADynamicsXY(); T)
    populate_vectors!(lsd2, N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy)

    gpu_csc, gpu_csr = if HAS_CUDA
        ins_cu   = gpu_inputs(ins)
        N_cu, N_ab_cu, ux_cu, uy_cu, taud_acx_cu, taud_acy_cu, β_acx_cu, β_acy_cu = ins_cu

        lsd2_csc = gpu_solver(lsd2)
        populate_vectors!(lsd2_csc, N_cu, N_ab_cu, ux_cu, uy_cu,
                          taud_acx_cu, taud_acy_cu, β_acx_cu, β_acy_cu)

        lsd2_csr = gpu_solver_csr(lsd2)
        populate_vectors!(lsd2_csr, N_cu, N_ab_cu, ux_cu, uy_cu,
                          taud_acx_cu, taud_acy_cu, β_acx_cu, β_acy_cu)
        CUDA.synchronize()

        (lsd2_csc, ins_cu), (lsd2_csr, ins_cu)
    else
        nothing, nothing
    end

    (; nx, ny, ins, lsd1, lsd2, gpu_csc, gpu_csr)
end

# ---------------------------------------------------------------------------
# Table 1: Assembly — populate_vectors!
# ---------------------------------------------------------------------------

w = 80
println("=" ^ w)
println(" Assembly — populate_vectors!")
println(" v1: COO fill + sparse(Ai,Aj,Av) | v2-CPU: KA kernel | v2-GPU: KA CUDABackend")
println("=" ^ w)
@printf("  %-18s %11s %11s %11s %9s %9s\n",
    "Grid (2×n DOF)", "v1 [μs]", "v2-CPU [μs]", "v2-GPU [μs]", "v2/v1", "GPU/CPU")
println("-" ^ w)

for (; nx, ny, ins, lsd1, lsd2, gpu_csc) in all_data
    N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy = ins

    t_fill  = (@b populate_vectors!(
        $lsd1, $N, $N_ab, $ux, $uy, $taud_acx, $taud_acy, $β_acx, $β_acy, DIVADynamicsXY()
    )).time
    t_assem = (@b sparse($lsd1.Ai, $lsd1.Aj, $lsd1.Av)).time
    t_v1    = t_fill + t_assem

    t_v2_cpu = (@b populate_vectors!(
        $lsd2, $N, $N_ab, $ux, $uy, $taud_acx, $taud_acy, $β_acx, $β_acy
    )).time

    if HAS_CUDA
        lsd2_gpu, ins_cu = gpu_csc
        N_cu, N_ab_cu, ux_cu, uy_cu, taud_acx_cu, taud_acy_cu, β_acx_cu, β_acy_cu = ins_cu

        t_v2_gpu = (@b begin
            populate_vectors!(
                $lsd2_gpu, $N_cu, $N_ab_cu, $ux_cu, $uy_cu,
                $taud_acx_cu, $taud_acy_cu, $β_acx_cu, $β_acy_cu,
            )
            CUDA.synchronize()
        end).time

        @printf("  %-18s %11.1f %11.1f %11.1f %8.1fx %8.1fx\n",
            "$(nx)×$(ny) ($(2*nx*ny))",
            t_v1 * 1e6, t_v2_cpu * 1e6, t_v2_gpu * 1e6,
            t_v1 / t_v2_cpu, t_v2_cpu / t_v2_gpu)
    else
        @printf("  %-18s %11.1f %11.1f %11s %8.1fx %9s\n",
            "$(nx)×$(ny) ($(2*nx*ny))",
            t_v1 * 1e6, t_v2_cpu * 1e6, "N/A",
            t_v1 / t_v2_cpu, "N/A")
    end
end

println("=" ^ w)

# ---------------------------------------------------------------------------
# Table 2: Solve — velocity!(lsd)
#
# v1: LinearProblem rebuilds sparse(Ai,Aj,Av) on every call, then solves.
# v2: LinearProblem uses the pre-built CSC matrix, then solves.
# GPU sparse direct solvers require cuSolver / additional LinearSolve extensions
# and are not benchmarked here.
# ---------------------------------------------------------------------------

println()
println("=" ^ w)
println(" Solve — velocity!(lsd)")
println(" v1: sparse(Ai,Aj,Av) + factorization | v2-CPU: pre-built matrix + factorization")
println(" v2-GPU: CUDSS direct LU via PagosCUDSSExt (CSR format)")
println("=" ^ w)
@printf("  %-18s %11s %11s %11s %8s %8s\n",
    "Grid (2×n DOF)", "v1 [ms]", "v2-CPU [ms]", "v2-GPU [ms]", "v2/v1", "GPU/CPU")
println("-" ^ w)

for (; nx, ny, ins, lsd1, lsd2, gpu_csr) in all_data
    N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy = ins

    # Repopulate so all solvers start from the same matrix/RHS.
    populate_vectors!(lsd1, N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy, DIVADynamicsXY())
    populate_vectors!(lsd2, N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy)

    t_v1 = (@b velocity!($lsd1)).time
    t_v2 = (@b velocity!($lsd2)).time

    if HAS_CUDA
        lsd2_gpu, _ = gpu_csr

        # Warmup to trigger JIT and lazy CUDSS factorization initialization.
        velocity!(lsd2_gpu)
        CUDA.synchronize()

        t_v2_gpu = (@b begin
            velocity!($lsd2_gpu)
            CUDA.synchronize()
        end).time

        @printf("  %-18s %11.3f %11.3f %11.3f %7.2fx %7.2fx\n",
            "$(nx)×$(ny) ($(2*nx*ny))",
            t_v1 * 1e3, t_v2 * 1e3, t_v2_gpu * 1e3,
            t_v1 / t_v2, t_v2 / t_v2_gpu)
    else
        @printf("  %-18s %11.3f %11.3f %11s %7.2fx %8s\n",
            "$(nx)×$(ny) ($(2*nx*ny))",
            t_v1 * 1e3, t_v2 * 1e3, "N/A",
            t_v1 / t_v2, "N/A")
    end
end

println("=" ^ w)
println()
println("Notes:")
println("  v1 assem       COO fill → sparse(Ai,Aj,Av) on every populate_vectors! call")
println("  v2-CPU assem   KA CPU() kernel → direct nzval write into pre-built CSC matrix")
println("  v2-GPU assem   KA CUDABackend() kernel, CuSparseMatrixCSC; includes CUDA.synchronize()")
println("  v1 solve       sparse() construction + UMFPACK symbolic+numeric factorization + ldiv!")
println("  v2-CPU solve   lu!/ldiv! with cached symbolic factorization (numeric only after 1st call)")
println("  v2-GPU solve   CUDSS numeric factorization only after 1st call (CSR, cached solver)")
HAS_CUDA || println("  GPU benchmarks skipped: no CUDA-capable device found (CUDA.functional() == false)")

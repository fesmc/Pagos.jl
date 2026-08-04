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

function gpu_solver(lsd2::LinearMomentumSolver2D)
    return LinearMomentumSolver2D(
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
function gpu_solver_csr(lsd2::LinearMomentumSolver2D)
    A_cpu = lsd2.A                              # SparseMatrixCSC on CPU
    Ai_csc, Aj_csc, _ = findnz(A_cpu)          # (row, col) in CSC nzval order
    csc_to_csr = invperm(sortperm(collect(zip(Ai_csc, Aj_csc))))
    perm_csr   = csc_to_csr[Array(lsd2.perm)]  # remap perm to CSR nzval positions
    return LinearMomentumSolver2D(
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

    lsd1 = LegacyLinearMomentumSolver2D(rp; T)
    populate_vectors!(lsd1, N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy, DIVADynamics())

    lsd2 = LinearMomentumSolver2D(rp, DIVADynamics(); T)
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
        $lsd1, $N, $N_ab, $ux, $uy, $taud_acx, $taud_acy, $β_acx, $β_acy, DIVADynamics()
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
# v1: rebuilds sparse(Ai,Aj,Av) on every call, then factorizes and solves via `\`.
# v2: reuses the pre-built CSC matrix, then factorizes and solves via lu!/ldiv!.
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
    populate_vectors!(lsd1, N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy, DIVADynamics())
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

# ---------------------------------------------------------------------------
# Pseudo-transient helpers
#
# pt_inputs   — build a NamedTuple of preallocated arrays (Matrix or CuArray)
# pt_one_iter! — one full PT velocity update; broadcast / circshift throughout
#                so the same code runs on CPU (Matrix) and GPU (CuArray)
# ---------------------------------------------------------------------------

function pt_inputs(nx, ny; H0 = 1000.0, μ0 = 1e5, β0 = 1e4, α = 1e-3,
                   ρ = 910.0, g = 9.81, dx = 5e3, T = Float64)
    H    = fill(T(H0), nx, ny)
    z_b  = T(-α) .* (reshape(T.(0:nx-1) .- T(nx ÷ 2), nx, 1) .* T(dx)) .* ones(T, 1, ny)
    mu   = fill(T(μ0), nx, ny)
    beta = fill(T(β0), nx, ny)
    z    = zeros(T, nx, ny)
    return (;
        H, z_b, mu, beta,
        beta_acx = copy(z), beta_acy = copy(z),
        N_ab     = copy(z),
        ux = copy(z), uy = copy(z),
        ux_old = copy(z), uy_old = copy(z),
        ux_x = copy(z), ux_y = copy(z),
        uy_x = copy(z), uy_y = copy(z),
        strainrate_xx = copy(z), strainrate_xy = copy(z), strainrate_yy = copy(z),
        shearstress_x = copy(z), shearstress_y = copy(z),
        basalstress_x = copy(z), basalstress_y = copy(z),
        drivingstress_x = copy(z), drivingstress_y = copy(z),
        dotvel_x = copy(z), dotvel_y = copy(z),
        prealloc = copy(z),
        dx = T(dx), dy = T(dx), rho_ice = T(ρ), g = T(g),
    )
end

function to_cu(s)
    return map(x -> x isa AbstractArray ? CuArray(x) : x, s)
end

function pt_one_iter!(s, dtau, theta_v)
    (; H, z_b, mu, beta, beta_acx, beta_acy, N_ab,
       ux, uy, ux_old, uy_old,
       ux_x, ux_y, uy_x, uy_y,
       strainrate_xx, strainrate_xy, strainrate_yy,
       shearstress_x, shearstress_y,
       basalstress_x, basalstress_y,
       drivingstress_x, drivingstress_y,
       dotvel_x, dotvel_y, prealloc,
       dx, dy, rho_ice, g) = s

    beta_acx .= 0.5 .* (beta .+ circshift(beta, (-1, 0)))
    beta_acy .= 0.5 .* (beta .+ circshift(beta, (0, -1)))
    @. prealloc = H * mu
    N_ab .= 0.25 .* (prealloc .+ circshift(prealloc, (-1, 0)) .+
                     circshift(prealloc, (0, -1)) .+ circshift(prealloc, (-1, -1)))

    velocitygradients!(ux_x, ux_y, uy_x, uy_y, ux, uy, dx, dy)
    scaledstrainrate!(strainrate_xx, strainrate_xy, strainrate_yy,
                      ux_x, ux_y, uy_x, uy_y, N_ab)
    shearstress!(shearstress_x, shearstress_y,
                 strainrate_xx, strainrate_xy, strainrate_yy, prealloc, dx, dy)

    ux_old .= ux
    uy_old .= uy
    basalstress!(basalstress_x, basalstress_y, beta_acx, beta_acy, ux_old, uy_old)
    drivingstress!(drivingstress_x, drivingstress_y, prealloc, rho_ice, g, H, z_b, dx, dy)

    Z = zero(eltype(H))
    @. dotvel_x = ifelse(H > 0, (shearstress_x - basalstress_x - drivingstress_x) / (rho_ice * H), Z)
    @. dotvel_y = ifelse(H > 0, (shearstress_y - basalstress_y - drivingstress_y) / (rho_ice * H), Z)

    pseudo_vel!(ux, ux_old, dotvel_x, dtau, theta_v)
    pseudo_vel!(uy, uy_old, dotvel_y, dtau, theta_v)
    return nothing
end

# Reset velocities and run full PT loop — repeated benchmark calls each start from zero.
function run_pt!(icesheet)
    icesheet.state.ux .= 0
    icesheet.state.uy .= 0
    pseudo_transient!(icesheet)
    return nothing
end

# ---------------------------------------------------------------------------
# Table 3: Single PT iteration — pt_one_iter! (CPU Matrix vs GPU CuArray)
# ---------------------------------------------------------------------------

const PT_SIZES = [(50, 50), (128, 128), (256, 256), (512, 512)]
const PT_DTAU    = T(1e-5)
const PT_THETA_V = T(0.6)

pt_all_data = map(PT_SIZES) do (nx, ny)
    s_cpu = pt_inputs(nx, ny; T)
    s_gpu = HAS_CUDA ? to_cu(s_cpu) : nothing
    (; nx, ny, s_cpu, s_gpu)
end

println()
println("=" ^ w)
println(" PT single iteration — pt_one_iter! (broadcast + circshift, GPU-compatible)")
println(" stagger_beta & vintegrated_viscosity via circshift; dotvel via ifelse broadcast")
println("=" ^ w)
@printf("  %-18s %11s %11s %8s\n", "Grid", "CPU [ms]", "GPU [ms]", "CPU/GPU")
println("-" ^ w)

for (; nx, ny, s_cpu, s_gpu) in pt_all_data
    t_cpu = (@b pt_one_iter!($s_cpu, $PT_DTAU, $PT_THETA_V)).time

    if HAS_CUDA
        pt_one_iter!(s_gpu, PT_DTAU, PT_THETA_V)   # warmup
        CUDA.synchronize()
        t_gpu = (@b begin
            pt_one_iter!($s_gpu, $PT_DTAU, $PT_THETA_V)
            CUDA.synchronize()
        end).time
        @printf("  %-18s %11.3f %11.3f %7.2fx\n",
            "$(nx)×$(ny)", t_cpu * 1e3, t_gpu * 1e3, t_cpu / t_gpu)
    else
        @printf("  %-18s %11.3f %11s %8s\n",
            "$(nx)×$(ny)", t_cpu * 1e3, "N/A", "N/A")
    end
end

println("=" ^ w)

# ---------------------------------------------------------------------------
# Table 4: Full PT convergence loop — pseudo_transient! (CPU, IceSheet path)
#          GPU path requires State to use AbstractMatrix; tracked as future work.
#          Velocity reset to zero before each sample so convergence is not trivial.
# ---------------------------------------------------------------------------

println()
println("=" ^ w)
println(" Full PT convergence — pseudo_transient! + velocity reset (CPU, IceSheet)")
println("=" ^ w)
@printf("  %-18s %11s\n", "Grid (2×n DOF)", "Time [ms]")
println("-" ^ w)

for (nx, ny) in PT_SIZES
    dx_pt = T(5e3)
    lx_pt = T(nx - 1) * dx_pt
    ly_pt = T(ny - 1) * dx_pt

    domain_pt  = Domain(T, lx_pt, ly_pt, dx_pt, dx_pt)
    state_pt   = State(domain_pt)
    params_pt  = Params{T}()
    options_pt = Options{T}(maxiter = 200, abstol = T(1e-8), printout_every = 99999)

    state_pt.H    .= T(1000.0)
    state_pt.z_b  .= T(-1e-3) .* domain_pt.X
    state_pt.mu   .= T(1e5)
    state_pt.beta .= T(1e4)

    icesheet_pt = IceSheet(state_pt, domain_pt, params_pt, options_pt)

    t = (@b run_pt!($icesheet_pt)).time
    @printf("  %-18s %11.3f\n", "$(nx)×$(ny) ($(2*nx*ny))", t * 1e3)
end

println("=" ^ w)
println()
println("Notes (PT):")
println("  pt_one_iter!       broadcast+circshift throughout; same code for CPU (Matrix) and GPU (CuArray)")
println("  run_pt!            resets ux=uy=0 then runs pseudo_transient! to convergence; β0=1e4 → ~10 iters")
println("  GPU full loop      not benchmarked: State{T} uses Matrix{T}; refactor to AbstractMatrix to enable")
HAS_CUDA || println("  GPU benchmarks skipped: no CUDA-capable device found (CUDA.functional() == false)")

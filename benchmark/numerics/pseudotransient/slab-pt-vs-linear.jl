# ---------------------------------------------------------------------------
# GPU benchmark: PseudoTransientSolver vs LinearMomentumSolver2D (CUDSS LU)
# on the uniform DIVA slab (same setup as test/mechanics/slab.jl).
#
# Timed quantities (per momentum solve, i.e. what one nonlinear/Picard
# iteration would pay):
#  - PT:      velocity!(mech, cst) from zero initial velocity, until the
#             max-norm velocity change < abstol.
#  - linear:  populate_vectors! (assembly) + velocity! (CUDSS numeric
#             factorization + triangular solves); the symbolic analysis is
#             done once outside the timed region, as in production use.
#
# CAVEAT: the uniform slab flatters PT. With uniform fields and strong basal
# friction the PT iteration count is small (printed below) and independent of
# grid size, while the direct factorization scales superlinearly. On problems
# with dominant membrane stresses the PT iteration count grows with resolution,
# so do not read the ratio below as a general statement — rerun on a realistic
# configuration before drawing conclusions.
#
# Run from the repo root:
#   JULIA_LOAD_PATH="$PWD:$PWD/benchmark:@stdlib" julia --startup-file=no \
#       benchmark/numerics/pseudotransient/slab-pt-vs-linear.jl
# ---------------------------------------------------------------------------

using Pagos, CUDA, CUDA.CUSPARSE, CUDSS, KernelAbstractions, SparseArrays, Printf

const T    = Float64
const SLAB = (H0 = 1000.0, μ0 = 1e5, β0 = 1e4, α = 1e-3)
const RHO, GRAV = 910.0, 9.81
const UB = RHO * GRAV * SLAB.H0 * SLAB.α / SLAB.β0    # analytical velocity
const DX = T(5e3)
const SIZES   = [(64, 64), (128, 128), (256, 256), (512, 512), (1024, 1024)]
const SAMPLES = 5                                     # min-of-N (display GPU is noisy)

# --- Pseudo-transient setup -------------------------------------------------

function pt_setup(nx, ny)
    grid   = RegularGrid(CUDABackend(), T, (nx - 1) * DX, (ny - 1) * DX, DX, DX)
    state  = MechanicState(grid)
    solver = PseudoTransientSolver(grid; maxiter = 10_000, abstol = 1e-8, ncheck = 10)
    mech   = Mechanics(state, grid, nothing, DIVAMomentumBalance(), solver,
                       nothing, nothing)
    state.topography.thickness .= T(SLAB.H0)
    state.topography.surface   .= T(SLAB.H0) .- T(SLAB.α) .* grid.x
    state.material.viscosity_depthaveraged .= T(SLAB.μ0)
    state.friction.beta_eff    .= T(SLAB.β0)
    return mech
end

function bench_pt(nx, ny)
    mech = pt_setup(nx, ny)
    cst  = Constants{T}()
    res  = velocity!(mech, cst)                      # warmup / compile
    t = Inf
    for _ in 1:SAMPLES
        fill!(mech.state.velocity.x, T(0))
        fill!(mech.state.velocity.y, T(0))
        t = min(t, @elapsed (res = velocity!(mech, cst)))
    end
    ux = Array(mech.state.velocity.x)
    @assert res.converged
    @assert all(x -> isapprox(x, UB; rtol = 1e-4), ux)
    return t, res.iterations
end

# --- Linear (CUDSS) setup ---------------------------------------------------

# The public constructor assembles the sparsity pattern on the CPU; CUDSS needs
# the matrix as CuSparseMatrixCSR, so remap perm from CSC to CSR nzval order
# and move all live arrays to the device.
function gpu_solver_csr(lsd)
    A_cpu      = lsd.A                                # SparseMatrixCSC on CPU
    Ai, Aj, _  = findnz(A_cpu)                        # (row, col) in CSC nzval order
    csc_to_csr = invperm(sortperm(collect(zip(Ai, Aj))))
    perm_csr   = csc_to_csr[Array(lsd.perm)]
    return LinearMomentumSolver2D(
        lsd.dynamics, lsd.nx, lsd.ny, lsd.dxdx_, lsd.dydy_, lsd.dxdy_,
        CuArray(lsd.u), CuArray(lsd.u0), CuArray(lsd.b),
        CuSparseMatrixCSR(A_cpu), CuArray(perm_csr),
        lsd.i_idx, lsd.j_idx, Ref{Any}(nothing),
    )
end

function lin_setup(nx, ny)
    grid = RegularGrid(T, (nx - 1) * DX, (ny - 1) * DX, DX, DX)
    lsd  = gpu_solver_csr(LinearMomentumSolver2D(grid, DIVAMomentumBalance()))
    f(v) = CUDA.fill(T(v), nx, ny)
    ins  = (f(SLAB.H0 * SLAB.μ0), f(SLAB.H0 * SLAB.μ0),          # N, N_ab
            f(0), f(0),                                          # ux, uy
            f(-RHO * GRAV * SLAB.H0 * SLAB.α), f(0),             # taud_acx, taud_acy
            f(SLAB.β0), f(SLAB.β0))                              # β_acx, β_acy
    return lsd, ins
end

function bench_lin(nx, ny)
    lsd, ins = lin_setup(nx, ny)
    populate_vectors!(lsd, ins...)                   # warmup: assembly kernels,
    velocity!(lsd)                                   # CUDSS analysis + factorization
    CUDA.synchronize()
    t = Inf
    for _ in 1:SAMPLES
        t = min(t, @elapsed begin
            populate_vectors!(lsd, ins...)
            velocity!(lsd)
            CUDA.synchronize()
        end)
    end
    ux = Array(lsd.u)[1:nx * ny]                     # ux block, n = (i-1)ny + j
    @assert all(x -> isapprox(x, UB; rtol = 1e-4), ux)
    return t
end

# --- Run ----------------------------------------------------------------------

@assert CUDA.functional()
println("GPU: ", CUDA.name(CUDA.device()), " | slab: ", SLAB, "\n")
w = 78
println("=" ^ w)
println(" DIVA uniform slab — PT (abstol 1e-8, ncheck 10) vs CUDSS LU, Float64")
println("=" ^ w)
@printf("  %-12s %9s %14s %14s %10s\n",
        "Grid", "PT iters", "PT [ms]", "linear [ms]", "lin/PT")
println("-" ^ w)
for (nx, ny) in SIZES
    t_pt, iters = bench_pt(nx, ny)
    t_lin       = bench_lin(nx, ny)
    @printf("  %-12s %9d %14.2f %14.2f %9.1fx\n",
            "$(nx)×$(ny)", iters, 1e3 * t_pt, 1e3 * t_lin, t_lin / t_pt)
end
println("=" ^ w)
println("\nlinear = assembly + numeric factorization + solve (symbolic analysis cached).")
println("See CAVEAT in the file header before generalizing the ratio.")

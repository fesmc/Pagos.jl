###############################################################
# Solvers
##############################################################

"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the dynamics solver via [`velocity`](@ref).

# Available subtypes:
"""
abstract type AbstractMomentumSolver end

# TODO
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via energy minimization.
"""
struct OptimMomentumSolver <: AbstractMomentumSolver
end

# TODO
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via an iterative linear solver (e.g., CG, GMRES).
"""
struct IterativeMomentumSolver <: AbstractMomentumSolver
end

"""
$(TYPEDSIGNATURES)

TODO: Solve the ice dynamics via wavelet methods.
"""
struct WaveletMomentumSolver <: AbstractMomentumSolver
end

# TODO
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via convolutional neural network (IGM style).
"""
struct ConvolutionalMomentumSolver <: AbstractMomentumSolver
end

# TODO
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via a transient solver (e.g., explicit time-stepping).
"""
struct TransientMomentumSolver <: AbstractMomentumSolver
end

"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via a pseudo-transient (PT) solver following Sandip et al. (2024).
The velocity field is relaxed in pseudo-time until the momentum balance is satisfied,
which only requires local (stencil) operations. All work arrays live on the backend of
the grid the solver is constructed from, so the same code runs on CPU and GPU.

Construct from a grid, overriding parameters selectively via keyword arguments:

```julia
solver = PseudoTransientSolver(grid)
solver = PseudoTransientSolver(grid; maxiter = 500, abstol = 1e-10)
```

# Fields:
 - `ndim1`, `ndim2`, `ndim3`: numerical-dimensionality constants of the PT time step
   (1D, 2D, 3D stencils).
 - `min_bulk_viscosity_ice`: lower bound of the bulk viscosity.
 - `muB`: bulk-to-shear viscosity ratio entering the PT time step.
 - `theta_v`: relaxation weight of the velocity update.
 - `theta_mu`: relaxation weight of the viscosity update.
 - `abstol`: convergence tolerance on the max-norm velocity change per iteration.
 - `maxiter`: maximum number of PT iterations.
 - `ncheck`: check convergence every `ncheck` iterations. The check is the only
   operation that forces the host to wait for the device (see
   [Asynchronous kernel launches](@ref)), so on GPU a larger value (10–50) keeps
   the launch pipeline busy at the cost of up to `ncheck - 1` extra iterations.
 - `printout_every`: print a convergence monitor every `printout_every` iterations
   (silent by default).
 - `dtau_scaling`: safety scaling of the PT time step.
 - `velocity_x_old`, `velocity_y_old`: previous velocity iterate.
 - `velocity_x_dt`, `velocity_y_dt`: pseudo-transient velocity rate.
"""
struct PseudoTransientSolver{T<:AbstractFloat, M} <: AbstractMomentumSolver
    ndim1::T
    ndim2::T
    ndim3::T
    min_bulk_viscosity_ice::T
    muB::T
    theta_v::T
    theta_mu::T
    abstol::T
    maxiter::Int
    ncheck::Int
    printout_every::Int
    dtau_scaling::T
    velocity_x_old::M
    velocity_y_old::M
    velocity_x_dt::M
    velocity_y_dt::M
end
Adapt.@adapt_structure PseudoTransientSolver

function PseudoTransientSolver(grid::RegularGrid;
    ndim1 = 2.1,
    ndim2 = 4.1,
    ndim3 = 6.1,
    min_bulk_viscosity_ice = 0.5,
    muB = 1e2,
    theta_v = 0.6,
    theta_mu = 0.1,
    abstol = 1e-8,
    maxiter = 100,
    ncheck = 1,
    printout_every = typemax(Int),
    dtau_scaling = 1,
)
    T = eltype(grid.x)
    backend = get_backend(grid.x)
    w() = KernelAbstractions.zeros(backend, T, grid.nx, grid.ny)
    return PseudoTransientSolver(
        T(ndim1), T(ndim2), T(ndim3), T(min_bulk_viscosity_ice), T(muB),
        T(theta_v), T(theta_mu), T(abstol), Int(maxiter), Int(ncheck),
        Int(printout_every), T(dtau_scaling), w(), w(), w(), w(),
    )
end

"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via a direct linear solver (e.g., sparse LU factorization).

# Improvements over `LegacyLinearMomentumSolver2D`:
 1. dynamics is a type parameter → dispatch on loop1!/loop2! is fully static, no need to thread a runtime `dynamics` argument through populate_vectors! layers.
 2. SparseMatrixCSC pre-allocated at construction (fixed sparsity pattern). The hot path writes directly to A.nzval via a precomputed COO→nzval index map, eliminating the sparse(Ai, Aj, Av) allocation on every solve.
 3. Single AI type parameter (i_idx and j_idx are always the same kind).
 4. VT/MT/PI type parameters for vector/matrix/perm arrays so the struct can hold GPU arrays (CuVector, CuSparseMatrix) without code changes. The populate_vectors! kernels are written with KernelAbstractions and run on whichever backend owns lsd.u.
"""
struct LinearMomentumSolver2D{
    DYN <: AbstractMomentumBalance,
    T   <: AbstractFloat,
    VT  <: AbstractVector,        # float vector type (u, u0, b)
    MT,                            # sparse matrix type (SparseMatrixCSC or CuSparseMatrix)
    PI  <: AbstractVector{Int},   # perm index vector type
    AI,
} <: AbstractMomentumSolver
    dynamics::DYN
    nx::Int
    ny::Int
    dxdx_::T
    dydy_::T
    dxdy_::T
    u::VT
    u0::VT
    b::VT
    A::MT
    perm::PI                      # COO fill order → A.nzval index
    i_idx::AI
    j_idx::AI
    solver_cache::Ref{Any}        # holds a cached direct solver (Nothing or backend-specific)
end

###############################################################
# Dispatch functions
###############################################################

# Functions to calculate velocity
function calc_F_integral(visc_eff, H_ice, f_ice, zeta_aa, n)
    # TODO: not yet implemented.
    error("calc_F_integral is not yet implemented")
end

# function velocity( solver::LinearSolver, dynamics::DIVA)
# end

function vertically_integrated_viscosity!(N, H, μ)
    N .= H .* μ
    return nothing
end

function stagger()
end

function fill_pattern!(Ai, Aj, nx, ny, i_idx, j_idx)
    k = 0
    for i in 1:nx, j in 1:ny
        im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
        nr = _ij2n_ux(i, j, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(ip1, j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(i,   j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(im1, j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(i,   jp1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(i,   jm1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(i,   j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(ip1, j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(ip1, jm1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(i,   jm1, nx, ny)
    end
    for i in 1:nx, j in 1:ny
        im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
        nr = _ij2n_uy(i, j, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(i,   jp1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(i,   j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(i,   jm1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(ip1, j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(im1, j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(i,   jp1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(i,   j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(im1, jp1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(im1, j,   nx, ny)
    end
    return nothing
end

function coo_to_nzval_idx(Ai, Aj, A::SparseMatrixCSC)
    perm = Vector{Int}(undef, length(Ai))
    for k in eachindex(Ai)
        c  = Aj[k]
        r  = Ai[k]
        lo = A.colptr[c]
        hi = A.colptr[c+1] - 1
        perm[k] = searchsortedfirst(A.rowval, r, lo, hi, Base.Order.Forward)
    end
    return perm
end

function LinearMomentumSolver2D(grid::RegularGrid, dynamics::DYN, backend = CPU()) where {DYN <: AbstractMomentumBalance}
    T      = eltype(grid.x)
    nx, ny = grid.nx, grid.ny
    dx, dy = T(grid.dx), T(grid.dy)
    dxdx_  = 1 / (dx * dx)
    dydy_  = 1 / (dy * dy)
    dxdy_  = 1 / (dx * dy)
    n_sprs = 18 * nx * ny  # 9 nonzeros/row × 2 equations (ux, uy)
    n_u    = 2 * nx * ny

    i_idx = PeriodicIndexing(1, nx)
    j_idx = PeriodicIndexing(1, ny)

    # Pattern and permutation are always computed on CPU (one-time cost).
    Ai_cpu   = zeros(Int, n_sprs)
    Aj_cpu   = zeros(Int, n_sprs)
    fill_pattern!(Ai_cpu, Aj_cpu, nx, ny, i_idx, j_idx)
    A_cpu    = sparse(Ai_cpu, Aj_cpu, ones(T, n_sprs), n_u, n_u)
    perm_cpu = coo_to_nzval_idx(Ai_cpu, Aj_cpu, A_cpu)

    # Allocate live arrays on the target backend.
    # For GPU (e.g. CUDABackend()), KernelAbstractions.zeros returns CuVector and
    # the sparse matrix should be adapted via CUDA.CUSPARSE.CuSparseMatrixCSC(A_cpu).
    u    = KernelAbstractions.zeros(backend, T,   n_u)
    u0   = KernelAbstractions.zeros(backend, T,   n_u)
    b    = KernelAbstractions.zeros(backend, T,   n_u)
    perm = KernelAbstractions.zeros(backend, Int, n_sprs)
    perm .= perm_cpu   # works for both CPU (no-op copy) and GPU (H→D transfer)

    # A stays as SparseMatrixCSC for the CPU default; for GPU, adapt before passing
    # to the inner constructor, e.g.:
    #   A = CUDA.CUSPARSE.CuSparseMatrixCSC(A_cpu)
    return LinearMomentumSolver2D(dynamics, nx, ny, dxdx_, dydy_, dxdy_, u, u0, b, A_cpu, perm, i_idx, j_idx, Ref{Any}(nothing))
end

@kernel function _assemble_ux!(nzval, perm, u0, b, nx, ny,
                               N, N_ab, ux, taud_acx, β_acx, β_acy,
                               dxdx_, dydy_, dxdy_, i_idx, j_idx, dynamics)
    i, j = @index(Global, NTuple)
    im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
    nr = _ij2n_ux(i, j, nx, ny)
    @inbounds u0[nr] = ux[i, j]
    @inbounds b[nr]  = taud_acx[i, j]
    v = loop1_coeffs(im1, i, ip1, jm1, j, jp1, dxdx_, dydy_, dxdy_,
                     N, N_ab, β_acx, β_acy, dynamics)
    k0 = 9 * ((i - 1) * ny + (j - 1))
    @inbounds for s in 1:9
        nzval[perm[k0 + s]] = v[s]
    end
end

@kernel function _assemble_uy!(nzval, perm, u0, b, nx, ny,
                               N, N_ab, uy, taud_acy, β_acx, β_acy,
                               dxdx_, dydy_, dxdy_, i_idx, j_idx, dynamics)
    i, j = @index(Global, NTuple)
    im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
    nr = _ij2n_uy(i, j, nx, ny)
    @inbounds u0[nr] = uy[i, j]
    @inbounds b[nr]  = taud_acy[i, j]
    v = loop2_coeffs(im1, i, ip1, jm1, j, jp1, dxdx_, dydy_, dxdy_,
                     N, N_ab, β_acx, β_acy, dynamics)
    k0 = 9 * nx * ny + 9 * ((i - 1) * ny + (j - 1))
    @inbounds for s in 1:9
        nzval[perm[k0 + s]] = v[s]
    end
end

"""
$(TYPEDSIGNATURES)

Assemble the linear system of the [`LinearMomentumSolver2D`](@ref) from the current dynamic
state `dyn_now`: fill the sparse matrix `A` (viscosity and basal-drag coefficients) and the
right-hand side `b` (driving stress), and seed the initial guess `u0` from the current
velocity. Mutates the solver buffers in place; the low-level method takes the unpacked fields
directly.
"""
function populate_vectors!(lsd::LinearMomentumSolver2D, dyn_now)
    (; N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy) = dyn_now
    populate_vectors!(lsd, N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy)
    return nothing
end

function populate_vectors!(
    lsd::LinearMomentumSolver2D,
    N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy,
)
    (; A, perm, u0, b, i_idx, j_idx, nx, ny, dxdx_, dydy_, dxdy_, dynamics) = lsd
    backend  = get_backend(lsd.u)
    kern_ux  = _assemble_ux!(backend)
    kern_uy  = _assemble_uy!(backend)

    nzval = nonzeros(A)
    kern_ux(
        nzval, perm, u0, b, nx, ny,
        N, N_ab, ux, taud_acx, β_acx, β_acy,
        dxdx_, dydy_, dxdy_, i_idx, j_idx, dynamics;
        ndrange = (nx, ny),
    )
    kern_uy(
        nzval, perm, u0, b, nx, ny,
        N, N_ab, uy, taud_acy, β_acx, β_acy,
        dxdx_, dydy_, dxdy_, i_idx, j_idx, dynamics;
        ndrange = (nx, ny),
    )
    return nothing
end

function LinearSolve.LinearProblem(lsd::LinearMomentumSolver2D)
    return LinearProblem(lsd.A, lsd.b; u0 = lsd.u)
end

function velocity!(lsd::LinearMomentumSolver2D)
    if lsd.solver_cache[] === nothing
        F = lu(lsd.A)
        lsd.solver_cache[] = F
        ldiv!(lsd.u, F, lsd.b)
    else
        _velocity_cached!(lsd.solver_cache[], lsd)
    end
    return nothing
end

# Function barrier: typed on F so lu! and ldiv! dispatch statically.
function _velocity_cached!(F, lsd::LinearMomentumSolver2D)
    lu!(F, lsd.A)
    ldiv!(lsd.u, F, lsd.b)
    return nothing
end

function velocity!(ux, uy, lsd::LinearMomentumSolver2D)
    u = lsd.u
    (; nx, ny) = lsd
    @inbounds for i in 1:nx, j in 1:ny
        ux[i, j] = u[_ij2n_ux(i, j, nx, ny)]
        uy[i, j] = u[_ij2n_uy(i, j, nx, ny)]
    end
    return nothing
end
"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the dynamics type via [`velocity`](@ref).

# Available subtypes:
"""
abstract type AbstractDynamics end

struct SIADynamics <: AbstractDynamics end
struct SSADynamics <: AbstractDynamics end
struct SIASSADynamics <: AbstractDynamics end
struct InertialSIASSADynamics <: AbstractDynamics end
struct DIVADynamics <: AbstractDynamics end
struct InertialDIVADynamics <: AbstractDynamics end
struct BlatterPattynDynamics <: AbstractDynamics end
struct StokesDynamics <: AbstractDynamics end

# Functions to calculate velocity 
function calc_F_integral(visc_eff,H_ice,f_ice,zeta_aa,n)
    # To do...
    return Fn
end

# function velocity( solver::LinearSolver, dynamics::DIVA)
# end

function vertically_integrated_viscosity!(N, H, μ)
    N .= H .* μ
    return nothing
end

function stagger()
end

struct ResolutionParameters{T}
    n::Int
    nx::Int
    ny::Int
    dx::T
    dy::T
    dxdx::T
    dydy::T
    dxdy::T
    dxdx_::T
    dydy_::T
    dxdy_::T
end

function ResolutionParameters(n, nx, ny, dx, dy; T=Float32)
    dx = T(dx)
    dy = T(dy)
    dxdx = dx * dx
    dydy = dy * dy
    dxdy = dx * dy
    dxdx_ = 1 / dxdx
    dydy_ = 1 / dydy
    dxdy_ = 1 / dxdy

    return ResolutionParameters(
        n, nx, ny,
        dx, dy,
        dxdx, dydy, dxdy,
        dxdx_, dydy_, dxdy_,
    )
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
        k+=1; Ai[k]=nr; Aj[k]=ij2n_uy(i,   j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=ij2n_uy(ip1, j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=ij2n_uy(ip1, jm1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=ij2n_uy(i,   jm1, nx, ny)
    end
    for i in 1:nx, j in 1:ny
        im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
        nr = ij2n_uy(i, j, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=ij2n_uy(i,   jp1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=ij2n_uy(i,   j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=ij2n_uy(i,   jm1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=ij2n_uy(ip1, j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=ij2n_uy(im1, j,   nx, ny)
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

function LinearDynamicsSolver2D(rp::ResolutionParameters, dynamics::DYN, backend = CPU(); T = Float32) where {DYN <: AbstractDynamics}
    (; n, nx, ny) = rp
    n_sprs = 2 * nx * ny * 3^n
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
    return LinearDynamicsSolver2D(dynamics, rp, u, u0, b, A_cpu, perm, i_idx, j_idx, Ref{Any}(nothing))
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
    nr = ij2n_uy(i, j, nx, ny)
    @inbounds u0[nr] = uy[i, j]
    @inbounds b[nr]  = taud_acy[i, j]
    v = loop2_coeffs(im1, i, ip1, jm1, j, jp1, dxdx_, dydy_, dxdy_,
                     N, N_ab, β_acx, β_acy, dynamics)
    k0 = 9 * nx * ny + 9 * ((i - 1) * ny + (j - 1))
    @inbounds for s in 1:9
        nzval[perm[k0 + s]] = v[s]
    end
end

function populate_vectors!(lsd::LinearDynamicsSolver2D, dyn_now)
    (; N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy) = dyn_now
    populate_vectors!(lsd, N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy)
    return nothing
end

function populate_vectors!(
    lsd::LinearDynamicsSolver2D,
    N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy,
)
    (; A, perm, u0, b, i_idx, j_idx, resolution_params, dynamics) = lsd
    (; nx, ny, dxdx_, dydy_, dxdy_) = resolution_params
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
    KernelAbstractions.synchronize(backend)
    return nothing
end

function LinearSolve.LinearProblem(lsd::LinearDynamicsSolver2D)
    return LinearProblem(lsd.A, lsd.b; u0 = lsd.u)
end

function velocity!(lsd::LinearDynamicsSolver2D)
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
function _velocity_cached!(F, lsd::LinearDynamicsSolver2D)
    lu!(F, lsd.A)
    ldiv!(lsd.u, F, lsd.b)
    return nothing
end

function velocity!(ux, uy, lsd::LinearDynamicsSolver2D)
    u = lsd.u
    (; nx, ny) = lsd.resolution_params
    @inbounds for i in 1:nx, j in 1:ny
        ux[i, j] = u[_ij2n_ux(i, j, nx, ny)]
        uy[i, j] = u[ij2n_uy(i, j, nx, ny)]
    end
    return nothing
end
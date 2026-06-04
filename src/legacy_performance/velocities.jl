
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via a direct linear solver (e.g., sparse LU factorization).
"""
struct LegacyLinearDynamicsSolver2D{
    RP,     # <: ResolutionParameters
    T,      # <: AbstractFloat
    AI1,    # <: AbstractIndexing
    AI2,    # <: AbstractIndexing
} <: AbstractDynamicsSolver
    resolution_params::RP
    n_terms::Int
    n_u::Int
    n_sprs::Int
    u::Vector{T}
    u0::Vector{T}
    b::Vector{T}
    Ai::Vector{Int}
    Aj::Vector{Int}
    Av::Vector{T}
    i_idx::AI1
    j_idx::AI2
end

function LegacyLinearDynamicsSolver2D(rp::ResolutionParameters; T=Float32)

    (; n, nx, ny) = rp
    n_terms = 3^n
    n_u     = 2*nx*ny
    n_sprs  = n_u*n_terms

    u = zeros(T, n_u)
    u0 = zeros(T, n_u)
    b = zeros(T, n_u)
    Ai = zeros(Int, n_sprs)
    Aj = zeros(Int, n_sprs)
    Av = zeros(T, n_sprs)
    i_idx = PeriodicIndexing(1, nx)
    j_idx = PeriodicIndexing(1, ny)

    return LegacyLinearDynamicsSolver2D(
        rp, n_terms, n_u, n_sprs,
        u, u0, b, Ai, Aj, Av,
        i_idx, j_idx,
    )
end


function populate_vectors!(lsd, dyn_now, dynamics)
    (; N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy) = dyn_now

    populate_vectors!(lsd, N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy, dynamics)
    return nothing
end

# TODO: RHS of Av should multiple dispatch on SSA and SIA too!

@inline function loop1_coeffs(::Any, i, ip1, jm1, j, ::Any, dxdx_, dydy_, dxdy_, N, N_ab, β_acx, ::Any, ::DIVADynamics)
    return (
         4dxdx_ * N[ip1, j],
        -4dxdx_ * (N[ip1, j] + N[i, j]) - dydy_ * (N_ab[i, j] + N_ab[i, jm1]) - β_acx[i, j],
         4dxdx_ * N[i, j],
         dydy_ * N_ab[i, j],
         dydy_ * N_ab[i, jm1],
        -2dxdy_ * N[i, j]   - dxdy_ * N_ab[i, j],
        -2dxdy_ * N[ip1, j] + dxdy_ * N_ab[i, j],
        -2dxdy_ * N[ip1, j] - dxdy_ * N_ab[i, jm1],
         2dxdy_ * N[i, j]   + dxdy_ * N_ab[i, jm1],
    )
end

@inline function loop2_coeffs(im1, i, ::Any, ::Any, j, jp1, dxdx_, dydy_, dxdy_, N, N_ab, ::Any, β_acy, ::DIVADynamics)
    return (
         4dydy_ * N[i, jp1],
        -4dydy_ * (N[i, jp1] + N[i, j]) - dxdx_ * (N_ab[i, j] + N_ab[im1, j]) - β_acy[i, j],
         4dydy_ * N[i, j],
         dxdx_ * N_ab[i, j],
         dxdx_ * N_ab[im1, j],
         2dxdy_ * N[i, jp1] + dxdy_ * N_ab[i, j],
        -2dxdy_ * N[i, j]   - dxdy_ * N_ab[i, j],
        -2dxdy_ * N[i, jp1] - dxdy_ * N_ab[im1, j],
         2dxdy_ * N[i, j]   + dxdy_ * N_ab[im1, j],
    )
end

function loop1!(v, im1, i, ip1, jm1, j, jp1, dxdx_, dydy_, dxdy_, N, N_ab, β_acx, β_acy, dyn::DIVADynamics)
    t = loop1_coeffs(im1, i, ip1, jm1, j, jp1, dxdx_, dydy_, dxdy_, N, N_ab, β_acx, β_acy, dyn)
    @inbounds for s in 1:9; v[s] = t[s]; end
    return nothing
end

function loop2!(v, im1, i, ip1, jm1, j, jp1, dxdx_, dydy_, dxdy_, N, N_ab, β_acx, β_acy, dyn::DIVADynamics)
    t = loop2_coeffs(im1, i, ip1, jm1, j, jp1, dxdx_, dydy_, dxdy_, N, N_ab, β_acx, β_acy, dyn)
    @inbounds for s in 1:9; v[s] = t[s]; end
    return nothing
end

function populate_vectors!(
    lsd::LegacyLinearDynamicsSolver2D,
    N,
    N_ab,
    ux,
    uy,
    taud_acx,
    taud_acy,
    β_acx,
    β_acy,
    dynamics,
)
    (; Ai, Aj, Av, u0, b, i_idx, j_idx, resolution_params) = lsd
    (; nx, ny, dxdx_, dydy_, dxdy_) = resolution_params

    populate_vectors!(
        Ai, Aj, Av, u0, b,
        nx, ny,
        N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy,
        dxdx_, dydy_, dxdy_,
        i_idx, j_idx,
        dynamics,
    )
    return nothing
end


function populate_vectors!(
    Ai, Aj, Av, u0, b,
    nx, ny,
    N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy,
    dxdx_, dydy_, dxdy_,
    i_idx, j_idx,
    dynamics,
)

    k = 0
    v = zeros(eltype(Av), 9)
    @inbounds for i in 1:nx, j in 1:ny

        im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)

        # Set the row in matrix A that the equation is being defined for:
        nr = _ij2n_ux(i, j, nx, ny)  
        u0[nr] = ux[i, j]
        b[nr] = taud_acx[i, j]
        loop1!(v, im1, i, ip1, jm1, j, jp1, dxdx_, dydy_, dxdy_, N, N_ab, β_acx, β_acy, dynamics::DIVADynamics)

        # -- vx terms --
        # vx(i+1, j)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_ux(ip1, j,nx,ny)
        Av[k] = v[1]
        

        # vx(i, j)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_ux(i, j,nx,ny)
        Av[k] = v[2]

        # vx(i-1, j)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_ux(im1, j,nx,ny)
        Av[k] = v[3]

        # vx(i, j+1)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_ux(i, jp1,nx,ny)
        Av[k] = v[4]

        # vx(i, j-1)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_ux(i, jm1,nx,ny)
        Av[k] = v[5]

        # -- vy terms -- 
        # vy(i, j)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_uy(i, j,nx,ny)
        Av[k] = v[6]

        # vy(i+1, j)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_uy(ip1, j,nx,ny)
        Av[k] = v[7]

        # vy(i+1, j-1)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_uy(ip1, jm1,nx,ny)
        Av[k] = v[8]
        
        # vy(i, j-1)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_uy(i, jm1,nx,ny)
        Av[k] = v[9]
        
    end

    for i in 1:nx, j in 1:ny
        im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)

        # Set the row in matrix A that the equation is being defined for:
        nr = _ij2n_uy(i, j, nx, ny)
        u0[nr] = uy[i, j]
        b[nr] = taud_acy[i, j]
        loop2!(v, im1, i, ip1, jm1, j, jp1, dxdx_, dydy_, dxdy_, N, N_ab, β_acx, β_acy, dynamics::DIVADynamics)

        # -- uy terms -- 
        # uy(i, j+1)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_uy(i, jp1,nx,ny)
        Av[k] = v[1]

        # uy(i, j)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_uy(i, j,nx,ny)
        Av[k] = v[2]

        # uy(i, j-1)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_uy(i, jm1,nx,ny)
        Av[k] = v[3]

        # uy(i+1, j)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_uy(ip1, j,nx,ny)
        Av[k] = v[4]

        # uy(i-1, j)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_uy(im1, j,nx,ny)
        Av[k] = v[5]

        # -- ux terms -- 
        # ux(i, j+1)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_ux(i, jp1,nx,ny)
        Av[k] = v[6]

        # ux(i, j)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_ux(i, j,nx,ny)
        Av[k] = v[7]

        # ux(i-1, j+1)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_ux(im1, jp1,nx,ny)
        Av[k] = v[8]
        
        # ux(i-1, j)
        k = k+1
        Ai[k] = nr
        Aj[k] = _ij2n_ux(im1, j,nx,ny)
        Av[k] = v[9]
    end
    return nothing
end


function velocity!(lsd::LegacyLinearDynamicsSolver2D; use_linsolve = true)
    if use_linsolve
        prob = LinearProblem(lsd)
        sol = solve(prob)
        lsd.u .= sol.u
    else
        lsd.u .= sparse(lsd.Ai, lsd.Aj, lsd.Av) \ lsd.b
    end
    return nothing
end

function velocity!(ux, uy, lsd::LegacyLinearDynamicsSolver2D)
    u = lsd.u
    (; nx, ny) = lsd.resolution_params
    @inbounds for i = 1:nx, j in 1:ny
        n1 = _ij2n_ux(i, j,nx,ny)
        ux[i, j] = u[n1]

        n2 = _ij2n_uy(i, j,nx,ny)
        uy[i, j] = u[n2]
    end
    return nothing
end


function LinearSolve.LinearProblem(lsd::LegacyLinearDynamicsSolver2D)
    return LinearProblem(sparse(lsd.Ai, lsd.Aj, lsd.Av), lsd.b; u0=lsd.u)
end
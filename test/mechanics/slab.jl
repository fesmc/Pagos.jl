using Pagos
using Test

# -----------------------------------------------------------------------
# Analytical ice-slab solution
#
# Uniform slab of thickness H0, effective viscosity μ0, basal friction β0,
# surface slope α, ice density ρ, and gravity g.
#
#   driving stress:       τd  = ρ g H0 α
#   DIVA (SSA) solution:  ub  = τd / β0          (basal / depth-averaged)
#   depth-averaged incl.
#   deformational flow:   ū   = ub + ρ g H0² α / (3 μ0)
# -----------------------------------------------------------------------

function slab_analytical(; H0, μ0, β0, α, ρ = 910.0, g = 9.81)
    τd = ρ * g * H0 * α
    ub = τd / β0
    ū  = ub + ρ * g * H0^2 * α / (3μ0)
    return (; τd, ub, ū)
end

# -----------------------------------------------------------------------
# Shared field setup for a uniform nx × ny slab.
#
# Periodic BCs are used so no boundary layers appear in the uniform field.
# N    = H0 μ0 — vertically integrated viscosity at i+½ faces (x-normal)
# N_ab = H0 μ0 — vertically integrated viscosity at j+½ faces (y-normal)
# Sign convention: b = ρgH ∂s/∂x < 0 for downslope flow in +x.
# -----------------------------------------------------------------------

function slab_fields(; H0, μ0, β0, α, ρ = 910.0, g = 9.81,
                      nx = 11, ny = 3, dx = 5e3, T = Float64)
    N     = fill(T(H0 * μ0),         nx, ny)
    N_ab  = fill(T(H0 * μ0),         nx, ny)
    β_acx = fill(T(β0),              nx, ny)
    β_acy = fill(T(β0),              nx, ny)
    τd_x  = fill(T(-ρ * g * H0 * α), nx, ny)
    τd_y  = fill(T(0),               nx, ny)
    ux    = zeros(T, nx, ny)
    uy    = zeros(T, nx, ny)
    grid  = RegularGrid(T, (nx - 1) * dx, (ny - 1) * dx, dx, dx)
    return N, N_ab, β_acx, β_acy, τd_x, τd_y, ux, uy, grid
end

# -----------------------------------------------------------------------
# Solvers
# -----------------------------------------------------------------------

function solve_slab_v1(; kw...)
    N, N_ab, β_acx, β_acy, τd_x, τd_y, ux, uy, grid = slab_fields(; kw...)
    rp  = Pagos.ResolutionParameters(2, grid.nx, grid.ny, grid.dx, grid.dy; T = eltype(ux))
    lsd = LegacyLinearMomentumSolver2D(rp; T = eltype(ux))
    populate_vectors!(lsd, N, N_ab, ux, uy, τd_x, τd_y, β_acx, β_acy, DIVAMomentumBalance())
    velocity!(lsd)
    velocity!(ux, uy, lsd)
    return ux, uy
end

function solve_slab_v2(; kw...)
    N, N_ab, β_acx, β_acy, τd_x, τd_y, ux, uy, grid = slab_fields(; kw...)
    lsd = LinearMomentumSolver2D(grid, DIVAMomentumBalance())
    populate_vectors!(lsd, N, N_ab, ux, uy, τd_x, τd_y, β_acx, β_acy)
    velocity!(lsd)
    velocity!(ux, uy, lsd)
    return ux, uy
end

# -----------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------

const SLAB_CASES = [
    (H0 = 1000.0, μ0 = 1e5,  β0 = 1e3,  α = 1e-3),
    (H0 =  500.0, μ0 = 4e5,  β0 = 30.0, α = 1e-3),
]

function check_slab(ux, uy, an)
    @test all(≈(an.ub, rtol = 1e-6), ux)
    @test all(≈(0.0,   atol = 1e-10 * abs(an.ub)), uy)
end

@testset "DIVA uniform slab — LegacyLinearMomentumSolver2D" begin
    for c in SLAB_CASES
        an = slab_analytical(; c...)
        ux, uy = solve_slab_v1(; c...)
        @testset "H0=$(c.H0) β0=$(c.β0)" begin
            check_slab(ux, uy, an)
        end
    end
end

@testset "DIVA uniform slab — LinearMomentumSolver2D" begin
    for c in SLAB_CASES
        an = slab_analytical(; c...)
        ux, uy = solve_slab_v2(; c...)
        @testset "H0=$(c.H0) β0=$(c.β0)" begin
            check_slab(ux, uy, an)
        end
    end
end

@testset "DIVA uniform slab — v2 repeated solves (lu! cached path)" begin
    # After the first velocity! call the solver_cache holds the LU object.
    # Subsequent calls reuse the symbolic factorization (lu!) and must still
    # converge to the correct solution when the matrix values change.
    c1, c2 = SLAB_CASES

    N1, N_ab1, β_acx1, β_acy1, τd_x1, τd_y1, ux1, uy1, grid = slab_fields(; c1...)
    N2, N_ab2, β_acx2, β_acy2, τd_x2, τd_y2, ux2, uy2, _    = slab_fields(; c2...)

    lsd = LinearMomentumSolver2D(grid, DIVAMomentumBalance())

    # First solve — exercises the lu() (cold) path and populates solver_cache.
    populate_vectors!(lsd, N1, N_ab1, ux1, uy1, τd_x1, τd_y1, β_acx1, β_acy1)
    velocity!(lsd)
    velocity!(ux1, uy1, lsd)
    an1 = slab_analytical(; c1...)
    @test lsd.solver_cache[] !== nothing
    @testset "1st call H0=$(c1.H0) β0=$(c1.β0)" begin
        check_slab(ux1, uy1, an1)
    end

    # Second solve — exercises the lu!() (cached symbolic) path with new field values.
    populate_vectors!(lsd, N2, N_ab2, ux2, uy2, τd_x2, τd_y2, β_acx2, β_acy2)
    velocity!(lsd)
    velocity!(ux2, uy2, lsd)
    an2 = slab_analytical(; c2...)
    @testset "2nd call H0=$(c2.H0) β0=$(c2.β0)" begin
        check_slab(ux2, uy2, an2)
    end
end

# -----------------------------------------------------------------------
# Operator properties on NON-uniform coefficients.
#
# Every test above solves with uy ≡ 0 on a uniform field, which is exactly
# the blind spot that let a sign error in loop1_coeffs' uy(i+1,j) entry
# (-2N[i+1,j]/(ΔxΔy) where the operator wants +2) survive: the term it
# corrupts multiplies uy, so uy ≡ 0 annihilates it, and a uniform N hides
# nothing else. These two checks need neither an analytic solution nor a
# reference implementation — they are properties the continuous SSA/DIVA
# operator has, so any correct discretization of it must have them too:
#
#   1. Self-adjointness: A is symmetric.
#   2. Rigid translation: with β = 0, a uniform velocity produces no force,
#      so [ux ≡ 1, uy ≡ 0] and [ux ≡ 0, uy ≡ 1] are both in the null space.
#
# Both need spatially varying N/N_ab to bite, hence the smooth fields below.
# -----------------------------------------------------------------------

using LinearAlgebra: norm, I

@testset "DIVA operator properties — non-uniform coefficients" begin
    nx, ny, dx = 12, 9, 5.0e3
    grid = RegularGrid(Float64, (nx - 1) * dx, (ny - 1) * dx, dx, dx)

    N     = [1e8 * (2 + sinpi(2i / nx) * cospi(2j / ny)) for i in 1:nx, j in 1:ny]
    N_ab  = [1e8 * (2 + cospi(2i / nx) * sinpi(2j / ny)) for i in 1:nx, j in 1:ny]
    β_acx = zeros(nx, ny)          # no drag: rigid translation must cost nothing
    β_acy = zeros(nx, ny)
    z     = zeros(nx, ny)

    lsd = LinearMomentumSolver2D(grid, DIVAMomentumBalance())
    populate_vectors!(lsd, N, N_ab, z, z, z, z, β_acx, β_acy)
    A = lsd.A
    n = nx * ny

    @test norm(A - A', Inf) ≤ 1e-12 * norm(A, Inf)

    for (name, u) in ("ux ≡ 1" => vcat(ones(n), zeros(n)),
                      "uy ≡ 1" => vcat(zeros(n), ones(n)))
        @testset "rigid translation $name" begin
            @test norm(A * u, Inf) ≤ 1e-12 * norm(A, Inf)
        end
    end
end

# The pseudo-transient solver used to be covered here too, against the same analytic slab
# on a collocated `RegularGrid`. That path is retired; its C-grid successor is tested in
# `pseudotransient_staggered.jl`, which carries every case this file had (uniform slab,
# ncheck > 1, Float32) plus the ones the collocated code never supported.

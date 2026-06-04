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
    rp    = ResolutionParameters(2, nx, ny, dx, dx; T)
    return N, N_ab, β_acx, β_acy, τd_x, τd_y, ux, uy, rp
end

# -----------------------------------------------------------------------
# Solvers
# -----------------------------------------------------------------------

function solve_slab_v1(; kw...)
    N, N_ab, β_acx, β_acy, τd_x, τd_y, ux, uy, rp = slab_fields(; kw...)
    lsd = LegacyLinearDynamicsSolver2D(rp; T = eltype(ux))
    populate_vectors!(lsd, N, N_ab, ux, uy, τd_x, τd_y, β_acx, β_acy, DIVADynamics())
    velocity!(lsd)
    velocity!(ux, uy, lsd)
    return ux, uy
end

function solve_slab_v2(; kw...)
    N, N_ab, β_acx, β_acy, τd_x, τd_y, ux, uy, rp = slab_fields(; kw...)
    lsd = LinearDynamicsSolver2D(rp, DIVADynamics(); T = eltype(ux))
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

@testset "DIVA uniform slab — LegacyLinearDynamicsSolver2D" begin
    for c in SLAB_CASES
        an = slab_analytical(; c...)
        ux, uy = solve_slab_v1(; c...)
        @testset "H0=$(c.H0) β0=$(c.β0)" begin
            check_slab(ux, uy, an)
        end
    end
end

@testset "DIVA uniform slab — LinearDynamicsSolver2D" begin
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

    N1, N_ab1, β_acx1, β_acy1, τd_x1, τd_y1, ux1, uy1, rp = slab_fields(; c1...)
    N2, N_ab2, β_acx2, β_acy2, τd_x2, τd_y2, ux2, uy2, _  = slab_fields(; c2...)

    lsd = LinearDynamicsSolver2D(rp, DIVADynamics(); T = Float64)

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
# Pseudo-transient tests
# -----------------------------------------------------------------------

@testset "Pseudo-transient pure functions" begin
    T   = Float64
    ρ   = T(910.0)
    dx  = T(5e3)
    muB = T(100.0)
    μ   = fill(T(1e5), 3, 3)

    expected_dt = ρ * dx^2 / (4 * (1 + muB) * 4.1 * T(1e5))
    @test pseudo_dt(ρ, dx, μ, muB) ≈ expected_dt

    v     = zeros(T, 2, 2)
    v_old = ones(T, 2, 2)
    dv    = fill(T(2.0), 2, 2)
    pseudo_vel!(v, v_old, dv, T(0.5), T(0.6))
    @test all(≈(T(1.6)), v)   # 1 + 0.6 * 2 * 0.5 = 1.6

    H       = fill(T(1000.0), 2, 2)
    shear   = fill(T(0.0),    2, 2)
    basal   = fill(T(1e3),    2, 2)
    driving = fill(T(-8927.0), 2, 2)
    dv2     = zeros(T, 2, 2)
    dotvel!(dv2, shear, basal, driving, ρ, H, 2, 2)
    @test all(≈((0.0 - 1e3 - (-8927.0)) / (ρ * T(1000.0))), dv2)
end

# Large β0 → lambda = θ·dtau·β/(ρH) ≈ 0.9 → converges in ~10 iterations.
const PT_SLAB = (H0 = 1000.0, μ0 = 1e5, β0 = 1e4, α = 1e-3)

@testset "DIVA uniform slab — pseudo_transient!" begin
    T   = Float64
    c   = PT_SLAB
    an  = slab_analytical(; c...)

    nx, ny, dx = 11, 3, T(5e3)
    lx = (nx - 1) * dx
    ly = (ny - 1) * dx

    domain  = Domain(T, lx, ly, dx, dx)
    state   = State(domain)
    params  = Params{T}()
    options = Options{T}(maxiter = 50, abstol = 1e-8, printout_every = 1000)

    state.H    .= T(c.H0)
    state.z_b  .= T(-c.α) .* domain.X   # linear slope: ∂(H+z_b)/∂x = -α
    state.mu   .= T(c.μ0)
    state.beta .= T(c.β0)

    icesheet = IceSheet(state, domain, params, options)
    pseudo_transient!(icesheet)

    @testset "H0=$(c.H0) β0=$(c.β0)" begin
        @test all(≈(an.ub, rtol = 1e-5), state.ux)
        @test all(≈(0.0,   atol = 1e-8 * abs(an.ub)), state.uy)
    end
end

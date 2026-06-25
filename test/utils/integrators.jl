using Pagos
using Test
using LinearAlgebra

# Linear decay problem:  dx/dt = -x,  x(0) = x0  =>  x(t) = x0 * exp(-t)
decay!(du, x, p, t) = (@. du = -x)

# Integrate `x0` from 0 to `T_end` with a freshly-constructed method.
function integrate(ctor, x0, T_end)
    x  = copy(x0)
    ig = Integrator(decay!, x, nothing, ctor(x))
    step!(ig, T_end)
    return ig.x
end

@testset "integrators" begin
    x0    = [1.0, -2.0, 3.0]
    T_end = 1.0
    exact = x0 .* exp(-T_end)

    fixed = (
        ("Euler",       x -> Euler(x; dt = 1e-3),       1e-3),
        ("RungeKutta4", x -> RungeKutta4(x; dt = 1e-2), 1e-8),
    )
    adaptive = (
        ("BogackiShampine32", x -> BogackiShampine32(x; dt = 0.1, atol = 1e-8, rtol = 1e-8)),
        ("Tsitouras54",       x -> Tsitouras54(x;       dt = 0.1, atol = 1e-10, rtol = 1e-10)),
    )

    # RKL2 needs a spectral-radius bound; for dx/dt = -x the Jacobian is -I, so ρ = 1.
    rkl2 = ("RKL2", x -> RKL2(x; dt = 0.1, spectral_radius = 1.0))

    @testset "type stability ($name)" for (name, ctor) in
            (fixed..., (adaptive[1][1], adaptive[1][2]), (adaptive[2][1], adaptive[2][2]), rkl2)
        ig = Integrator(decay!, copy(x0), nothing, ctor(x0))
        @inferred step!(ig, T_end)
    end

    @testset "accuracy ($name)" for (name, ctor, tol) in fixed
        @test integrate(ctor, x0, T_end) ≈ exact rtol = tol
    end

    @testset "accuracy ($name)" for (name, ctor) in adaptive
        @test integrate(ctor, x0, T_end) ≈ exact rtol = 1e-5
    end

    @testset "convergence order" begin
        # Halving the step should reduce the error by ~2^order.
        err(ctor) = maximum(abs.(integrate(ctor, x0, T_end) .- exact))

        e_euler_h  = err(x -> Euler(x; dt = 1e-2))
        e_euler_h2 = err(x -> Euler(x; dt = 5e-3))
        @test e_euler_h / e_euler_h2 ≈ 2 atol = 0.2          # order 1

        e_rk4_h  = err(x -> RungeKutta4(x; dt = 1e-1))
        e_rk4_h2 = err(x -> RungeKutta4(x; dt = 5e-2))
        @test e_rk4_h / e_rk4_h2 ≈ 16 atol = 3               # order 4
    end

    @testset "adaptive respects tolerance" begin
        # Tighter tolerance ⇒ smaller (or equal) error.
        loose = integrate(x -> Tsitouras54(x; dt = 0.2, atol = 1e-4, rtol = 1e-4), x0, T_end)
        tight = integrate(x -> Tsitouras54(x; dt = 0.2, atol = 1e-10, rtol = 1e-10), x0, T_end)
        @test maximum(abs.(tight .- exact)) ≤ maximum(abs.(loose .- exact))
    end

    @testset "RKL2 accuracy and order" begin
        # Second order on the (non-stiff) decay problem, ρ = 1.
        @test integrate(x -> RKL2(x; dt = 0.1, spectral_radius = 1.0), x0, T_end) ≈ exact rtol = 5e-3

        err(dt) = maximum(abs.(integrate(x -> RKL2(x; dt = dt, spectral_radius = 1.0), x0, T_end) .- exact))
        @test err(0.1) / err(0.05) ≈ 4 atol = 1.0           # order 2
    end

    @testset "RKL2 automatic spectral-radius estimate" begin
        # Diagonal operator f(u) = -λ.*u  ⇒  J = diag(-λ), spectral radius = max(λ),
        # with a well-separated dominant entry so the power method converges cleanly.
        λv = Float64[1, 3, 100, 7, 2, 40, 5, 9]
        diagop!(du, u, q, t) = (@. du = -λv * u)
        u0 = ones(length(λv))
        m  = RKL2(u0; dt = 1.0)                       # auto mode (spectral_radius = nothing)
        @test m.auto
        copyto!(m.y0, u0); diagop!(m.f0, m.y0, nothing, 0.0)
        ρ = estimate_spectral_radius!(diagop!, nothing, 0.0, m)
        @test ρ / m.safety ≈ maximum(λv) rtol = 0.05  # recovers max|λ|, not a sub-dominant one

        # Auto mode also stays type-stable and accurate end-to-end (decay: ρ = 1).
        ig = Integrator(decay!, copy(x0), nothing, RKL2(x0; dt = 0.1))
        @inferred step!(ig, T_end)
        @test ig.x ≈ exact rtol = 5e-3
    end

    @testset "RKL2 super-time-stepping stability" begin
        # Stiff scalar decay dx/dt = -λx. Super-time-stepping buys *stability*, not
        # accuracy, over a superstep that is large relative to the explicit limit.
        λ = 100.0
        stiff!(du, x, p, t) = (@. du = -λ * x)

        # A single explicit Euler step with h·λ = 100 ≫ 2 is violently unstable.
        @test abs(1.0 + 1.0 * (-λ * 1.0)) > 1

        # RKL2 absorbs the same superstep by raising the stage count: bounded & decaying.
        run(h) = (xs = [1.0];
                  step!(Integrator(stiff!, xs, nothing, RKL2(xs; dt = h, spectral_radius = λ)), 1.0);
                  xs[1])
        r_coarse = run(1.0)     # 1 superstep over [0,1]
        r_fine   = run(0.1)     # 10 supersteps
        @test isfinite(r_coarse)
        @test 0 ≤ r_coarse < 1                              # no blow-up, monotone decay
        @test r_fine < r_coarse                             # resolving the timescale ⇒ → e^{-λ} ≈ 0
    end

    @testset "SSPRK accuracy, order, and SSP property" begin
        @test ssp_coefficient(SSPRK33(x0; dt_fe = 0.1)) == 1
        @test ssp_coefficient(SSPRK43(x0; dt_fe = 0.1)) == 2

        # Both are third order on the smooth decay problem (step = C·dt_fe).
        integ_ssp(ctor, dt_fe) = (ig = Integrator(decay!, copy(x0), nothing, ctor(x0; dt_fe = dt_fe));
                                  step!(ig, T_end); ig.x)
        @test integ_ssp(SSPRK33, 0.05)  ≈ exact rtol = 1e-2
        @test integ_ssp(SSPRK43, 0.025) ≈ exact rtol = 1e-2

        err33(dt_fe) = maximum(abs.(integ_ssp(SSPRK33, dt_fe) .- exact))
        @test err33(0.1) / err33(0.05) ≈ 8 atol = 2.5          # order 3

        @inferred step!(Integrator(decay!, copy(x0), nothing, SSPRK33(x0; dt_fe = 0.1)), T_end)
        @inferred step!(Integrator(decay!, copy(x0), nothing, SSPRK43(x0; dt_fe = 0.1)), T_end)

        # SSP/positivity: at h = C·dt_fe the SSP method keeps a positive state ≥ 0, where
        # forward Euler at the same h undershoots below zero. (dx/dt = -x; FE positive
        # only for h ≤ 1, so h = 2 breaks it; SSPRK43 with C=2, dt_fe=1 steps h=2 safely.)
        decE!(du, u, q, t) = (@. du = -u)
        xs = [1.0]
        step!(Integrator(decE!, xs, nothing, SSPRK43(xs; dt_fe = 1.0)), 2.0)
        @test xs[1] ≥ 0
        @test (1.0 + 2.0 * (-1.0)) < 0                         # forward Euler at h = 2 goes negative
    end
end

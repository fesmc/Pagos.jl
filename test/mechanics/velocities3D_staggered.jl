using Pagos
using Test

include("../test_helpers/chmy.jl")

# Stage 3 of the DIVA/SSA split (`roadmaps/chmy.md`, Phase 3, decisions 13-15): the 3D
# velocity reconstruction — Robinson et al. (2022)'s second step, "Eq. (16) is integrated
# vertically to find the 3D velocity", run *after* `pseudo_transient!` converges, never
# inside it.

@testset "3D velocity reconstruction (C-grid staggered)" begin

    # Isolates the reconstruction kernel itself from F₂'s already-documented O(Δζ²)
    # quadrature error (Stage 1): using the code's *own* computed F₂ (not the pure
    # H/(3µ) closed form) to build the reference u_b makes this an exactness check on the
    # per-layer antiderivative g(σ) = σ - σ²/2 alone, at whatever nz — not a convergence
    # story. Uniform µ makes the per-layer piecewise-constant assumption exactly true, so
    # roundoff is the right bar.
    @testset "DIVA reconstruction: exact given the code's own F₁/F₂ (uniform µ)" begin
        H0, μ0, β0 = 1000.0, 1e5, 1e4
        ūx0, ūy0 = 20.0, -7.0
        lay  = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 12))
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, lay)
        rt   = Runtime(grid)
        mech = MechanicState(grid)

        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
        fill_analytic!(mech.friction.beta, rt.grid2d, (x, y) -> β0)
        fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d, (x, y) -> ūx0)
        fill_analytic!(mech.velocity.depthaverage_y, rt.grid2d, (x, y) -> ūy0)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0)

        solver = PseudoTransientSolver(grid; viscosity_continuation = NoViscosityContinuation())
        diva_update!(mech, solver, rt)
        velocities3D!(mech, rt, DIVAMomentumBalance())

        F2code = interior(mech.material.viscosity_integral_2)[4, 4, 1]
        F1code = interior(mech.material.viscosity_integral_1)[4, 4, 1]
        ub_x = ūx0 / (1 + β0 * F2code)
        ub_y = ūy0 / (1 + β0 * F2code)

        for k in 1:grid.nz
            ζ = zcenter(rt.grid, k)
            g(σ) = σ - σ^2 / 2
            F1p = (H0 / μ0) * g(ζ)   # antiderivative from 0, g(0) = 0
            @test all(≈(ub_x * (1 + β0 * F1p), rtol = 1e-12),
                     interior(mech.velocity.x)[:, :, k])
            @test all(≈(ub_y * (1 + β0 * F1p), rtol = 1e-12),
                     interior(mech.velocity.y)[:, :, k])
        end
        @test all(≈(ub_x * (1 + β0 * F1code), rtol = 1e-12), interior(mech.velocity.surface_x))
        @test all(≈(ub_y * (1 + β0 * F1code), rtol = 1e-12), interior(mech.velocity.surface_y))
    end

    # The full pipeline (diva_update! + reconstruction) against the *true* exact closed
    # form, for a genuinely depth-varying µ(σ) = µ0/(1+σ) (Stage 1's own test profile).
    # Unlike the test above, this does *not* subtract out F₂'s quadrature error — it is
    # the black-box "does this converge at the rate the underlying integrals promise"
    # check, and the promised rate is F₂'s second order (the per-layer piecewise-constant
    # µ discretization has that same order for a smoothly-varying profile).
    @testset "DIVA reconstruction: 2nd-order convergence vs. the true closed form" begin
        H0, μ0, β0 = 1000.0, 1e5, 1e4
        ūx0 = 20.0
        F2exact = 5H0 / (12μ0)     # Stage 1's _F2_VARYING, μ(σ) = μ0/(1+σ)
        ub_exact = ūx0 / (1 + β0 * F2exact)

        errs = Float64[]
        for nz in (8, 16, 32)
            lay  = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, nz))
            grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, lay)
            rt   = Runtime(grid)
            mech = MechanicState(grid)

            fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
            fill_analytic!(mech.friction.beta, rt.grid2d, (x, y) -> β0)
            fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d, (x, y) -> ūx0)
            fill_analytic!(mech.velocity.depthaverage_y, rt.grid2d, (x, y) -> 0.0)
            fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0 / (1 + ζ))

            solver = PseudoTransientSolver(grid; viscosity_continuation = NoViscosityContinuation())
            diva_update!(mech, solver, rt)
            velocities3D!(mech, rt, DIVAMomentumBalance())

            maxerr = 0.0
            for k in 1:grid.nz
                ζ = zcenter(rt.grid, k)
                # (1-σ)/μ(σ) = (1-σ)(1+σ)/μ0 = (1-σ²)/μ0 ⟹ ∫₀^ζ = ζ - ζ³/3.
                F1p_exact = (H0 / μ0) * (ζ - ζ^3 / 3)
                ux_exact = ub_exact * (1 + β0 * F1p_exact)
                got = interior(mech.velocity.x)[4, 4, k]
                maxerr = max(maxerr, abs(got - ux_exact) / abs(ux_exact))
            end
            push!(errs, maxerr)
        end
        @test all(>(1.7), convergence_rates(errs))
        @test errs[end] < 1e-3
    end

    # Structural cross-check: the reconstruction's own per-layer running sum (via the exact
    # antiderivative g) and `viscosity_integrals!`'s independently-coded F₁ (via the
    # midpoint-rule formula) are two different code paths computing the *same* quantity —
    # algebraically identical for a linear integrand, but not guaranteed bit-identical
    # given different floating-point evaluation order. `≈`, not `==`.
    @testset "surface velocity's F₁ matches viscosity_integrals!'s F₁ independently" begin
        lay  = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 10))
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, lay)
        rt   = Runtime(grid)
        mech = MechanicState(grid)
        H0, μ0, β0 = 1000.0, 1e5, 1e4

        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
        fill_analytic!(mech.friction.beta, rt.grid2d, (x, y) -> β0)
        fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d, (x, y) -> 15.0)
        fill_analytic!(mech.velocity.depthaverage_y, rt.grid2d, (x, y) -> 0.0)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> 1e5 / (1 + ζ))

        F1_ref = Field(rt.arch, rt.grid2d, (Center(), Center(), Center()))
        F2_ref = Field(rt.arch, rt.grid2d, (Center(), Center(), Center()))
        viscosity_integrals!(F1_ref, F2_ref, mech, rt)

        solver = PseudoTransientSolver(grid; viscosity_continuation = NoViscosityContinuation())
        diva_update!(mech, solver, rt)
        @test all(isapprox.(interior(mech.material.viscosity_integral_1), interior(F1_ref);
                           rtol = 1e-10))
    end

    # Physical sanity, independent of exact numbers: the DIVA velocity profile is a
    # monotone interpolation from `u_b` (frozen bed retards flow) to `u_s` (free surface
    # slides fastest) whenever β, µ > 0 — Eq. 16's `(1+βF₁(ζ))` factor is increasing in ζ.
    # A sign error anywhere in the reconstruction (e.g. the wrong side of `g`) would show
    # up here even if it happened to pass a narrower numeric check by coincidence.
    @testset "profile is monotone base → surface" begin
        lay  = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 10))
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, lay)
        rt   = Runtime(grid)
        mech = MechanicState(grid)
        H0, μ0, β0 = 1000.0, 1e5, 1e4

        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
        fill_analytic!(mech.friction.beta, rt.grid2d, (x, y) -> β0)
        fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d, (x, y) -> 20.0)
        fill_analytic!(mech.velocity.depthaverage_y, rt.grid2d, (x, y) -> 0.0)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0)

        solver = PseudoTransientSolver(grid; viscosity_continuation = NoViscosityContinuation())
        diva_update!(mech, solver, rt)
        velocities3D!(mech, rt, DIVAMomentumBalance())

        prof = interior(mech.velocity.x)[4, 4, :]
        @test issorted(prof)
        F2 = interior(mech.material.viscosity_integral_2)[4, 4, 1]
        ub = 20.0 / (1 + β0 * F2)
        @test ub < prof[1]
        @test prof[end] < interior(mech.velocity.surface_x)[4, 4, 1]
    end

    # Ice-free / masked columns get zero throughout, not a division artifact — the same
    # containment pattern `viscosity_integrals!` already established.
    @testset "H = 0 gives zero, not NaN/Inf" begin
        lay  = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 8))
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, lay)
        rt   = Runtime(grid)
        mech = MechanicState(grid)

        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> x < 0 ? 0.0 : 800.0)
        fill_analytic!(mech.friction.beta, rt.grid2d, (x, y) -> 1e4)
        fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d, (x, y) -> 10.0)
        fill_analytic!(mech.velocity.depthaverage_y, rt.grid2d, (x, y) -> 0.0)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> 1e5)

        solver = PseudoTransientSolver(grid; viscosity_continuation = NoViscosityContinuation())
        diva_update!(mech, solver, rt)
        velocities3D!(mech, rt, DIVAMomentumBalance())

        @test !any(isnan, interior(mech.velocity.x))
        @test !any(isinf, interior(mech.velocity.x))
        for i in axes(interior(mech.velocity.x), 1)
            x, _, _ = coord(rt.grid2d, location(mech.velocity.depthaverage_x), i, 1, 1)
            x < 0 && @test all(==(0.0), interior(mech.velocity.x)[i, :, :])
        end
    end

    # SSA is plug flow: every layer, and the surface, exactly equal the depth-averaged
    # solve — no shear, no viscosity integrals touched at all (decision 15).
    @testset "SSA reconstruction: exact plug flow" begin
        lay  = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 6))
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, lay)
        rt   = Runtime(grid)
        mech = MechanicState(grid)

        fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d, (x, y) -> 3.0x - 1.5y)
        fill_analytic!(mech.velocity.depthaverage_y, rt.grid2d, (x, y) -> 0.7x + 2.0y)

        velocities3D!(mech, rt, SSAMomentumBalance())

        ux2d = interior(mech.velocity.depthaverage_x)
        uy2d = interior(mech.velocity.depthaverage_y)
        for k in 1:grid.nz
            @test interior(mech.velocity.x)[:, :, k] == ux2d[:, :, 1]
            @test interior(mech.velocity.y)[:, :, k] == uy2d[:, :, 1]
        end
        @test interior(mech.velocity.surface_x) == ux2d
        @test interior(mech.velocity.surface_y) == uy2d
    end
end

using Pagos
using Test
using LinearAlgebra: Tridiagonal

include("../test_helpers/chmy.jl")

# Phase 2 of `roadmaps/blatter-pattyn.md`: vertical-implicit line relaxation. The claim is
# narrow and sharp — *the same fixed point as Phase 1, reached in an iteration count that does
# not grow with `nz`* — so the tests below are exactly those two statements plus the two
# invariants the implementation rests on (the Δτ bound loses its vertical term; the tridiagonal
# reduces to `pseudo_vel!` when there is no vertical operator to invert).
#
# Phase 1's own lesson (a) applies here in full: every uniform-slab check runs on geometry
# where `node_active` and `node_fully_active` agree everywhere, and the one Phase 1 bug unit
# tests missed lived entirely in the difference. So the fixed-point comparison below is run on
# a *masked* domain with lateral variation, not on a clean slab.

# Lateral variation + genuine vertical shear + a margin, small enough to solve twice quickly.
# `µ` is zero off-ice on purpose (the case the strict-mask rule was reached for), so this also
# exercises `_mu_acxz`'s one-sided interpolation through the tridiagonal assembly.
function sheared_margin_case(; nz, nx = 12, ny = 6, dx = 2e3)
    layering = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, nz))
    grid = StaggeredGrid(Float64, nx * dx, ny * dx, dx, dx, layering)
    rt   = Runtime(grid)
    mech = MechanicState(grid)
    topo = TopographicState(grid)

    m = nx - 3                       # ice out to i = m, ice-free beyond: one clean margin
    setdata!(topo.mask.is_ice, false)
    interior(topo.mask.is_ice)[1:m, :, 1] .= true
    mask = IceMask(topo.mask.is_ice)

    # Thickness and surface both vary in x *and* y, so the membrane terms are genuinely live
    # (a uniform slab kills them identically and would not test the horizontal operator at all).
    H(x, y) = 800.0 + 400.0 * sin(2π * x / (nx * dx)) + 150.0 * cos(2π * y / (ny * dx))
    s(x, y) = -1e-2 * x + 30.0 * sin(2π * y / (ny * dx))
    fill_analytic!(mech.topography.thickness, rt.grid2d, H)
    fill_analytic!(mech.topography.surface, rt.grid2d, s)
    # Strong basal-friction contrast across the domain: the bed row is where Phase 2's second
    # win lives, so the two treatments must agree across a β that spans four orders.
    fill_analytic!(mech.friction.beta_eff, rt.grid2d,
                   (x, y) -> 1e1 + 1e5 * exp(-((x - nx * dx / 3) / (2dx))^2))

    setdata!(mech.material.viscosity, 0.0)
    for k in 1:nz, j in -1:(ny + 2), i in 1:m
        # Depth-varying µ: a stiffer surface over softer basal ice, so `R₋ ≠ R₊` everywhere
        # and the tridiagonal's off-diagonals are genuinely asymmetric.
        mech.material.viscosity[i, j, k] = 5e7 * (1 + 2 * zcenter(rt.grid, k))
    end
    return (; grid, rt, mech, mask, nx, ny, nz, m)
end

function solve_case(case, vertical; abstol = 1e-7, maxiter = 60_000)
    (; grid, rt, mech, mask) = case
    momentum = BlatterPattynMomentumBalance()
    setdata!(mech.velocity.x, 0.0)
    setdata!(mech.velocity.y, 0.0)
    solver = PseudoTransientSolver(grid, momentum;
        maxiter, abstol, ncheck = 20,
        pseudo_timestep = GershgorinPseudoTimeStep(cfl = 0.99),
        convergence = ScaledResidual(),
        tuning = AutotunedDynamicRelaxation(),
        vertical_treatment = vertical)
    res = pseudo_transient!(mech, Constants{Float64}(), solver, rt, momentum, mask)
    return (; res, solver)
end

@testset "Blatter-Pattyn vertical-implicit line relaxation" begin
    cst = Constants{Float64}()
    layering(T = Float64; nz = 4) = CorrectedVerticalLayering(T, QuadraticSigmaTransform(T, nz))

    @testset "ImplicitVertical construction and guards" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, layering())
        vt = ImplicitVertical(grid)

        # Scratch is one column field per velocity component, at the component's own node.
        @test location(vt.thomas_x) === (Vertex(), Center(), Center())
        @test location(vt.thomas_y) === (Center(), Vertex(), Center())
        @test size(interior(vt.thomas_x)) == (grid.nx + 1, grid.ny, grid.nz)

        solver = PseudoTransientSolver(grid, BlatterPattynMomentumBalance();
                                       vertical_treatment = vt)
        @test solver.vertical_treatment === vt
        @test PseudoTransientSolver(grid, BlatterPattynMomentumBalance()
                                    ).vertical_treatment === ExplicitVertical()

        # The depth-averaged solver would never read it, so it says so rather than ignoring it.
        @test_throws ArgumentError PseudoTransientSolver(grid; vertical_treatment = vt)
        # ... and a strategy built from a different grid than the solver is a shape mismatch
        # that would otherwise only surface as a bounds error inside a kernel.
        other = StaggeredGrid(Float64, 16.0, 8.0, 1.0, 1.0, layering())
        @test_throws ArgumentError PseudoTransientSolver(
            grid, BlatterPattynMomentumBalance(); vertical_treatment = ImplicitVertical(other))
    end

    @testset "Δτ loses its vertical term, and only that" begin
        # The bound is the whole reason this phase exists: under `ImplicitVertical` the
        # vertical rows are inverted rather than stepped over, so `Δτ` must be set by
        # `Λ_horiz` alone. At Δx = 1 km against H = 100 m the vertical term dominates by
        # ~4 orders (the aspect ratio of `roadmaps/blatter-pattyn.md` §2), which is what makes
        # this measurable rather than a rounding difference.
        nx, ny, nz = 6, 4, 6
        H0, μ0, β0, dx = 100.0, 1e7, 1e3, 1e3
        grid = StaggeredGrid(Float64, nx * dx, ny * dx, dx, dx, layering(; nz))
        rt = Runtime(grid)
        mech = MechanicState(grid)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0)
        fill_analytic!(mech.friction.beta_eff, rt.grid2d, (x, y) -> β0)

        momentum = BlatterPattynMomentumBalance()
        exp_solver = PseudoTransientSolver(grid, momentum)
        imp_solver = PseudoTransientSolver(grid, momentum;
                                           vertical_treatment = ImplicitVertical(grid))
        pseudo_dt!(exp_solver, mech, cst, rt, momentum)
        pseudo_dt!(imp_solver, mech, cst, rt, momentum)

        @test all(>(0), interior(imp_solver.dtau_x))
        @test all(isfinite, interior(imp_solver.dtau_x))
        # Per row, not min-against-max: the gain *is* `Λ_total/Λ_horiz` at that row, and it
        # varies down the column with the sigma spacing (thin end layers gain the most).
        gain = interior(imp_solver.dtau_x) ./ interior(exp_solver.dtau_x)
        @test all(>(50), gain)
        @test maximum(gain) > 1e3

        # And it must be *exactly* the horizontal row sum, not merely a larger number: at
        # uniform µ with dx = dy the membrane rows give Λ = (16µ + 8µ + 4µ + 4µ)/ρ per §2.1
        # (P₋+P₊ = Q₋+Q₊ = 2µ), i.e. Δτ = 2·cfl/Λ.
        Λh = (8 * 2μ0 / dx^2 + 4 * 2μ0 / dx^2 + 2 * 2μ0 / dx^2 + 2 * 2μ0 / dx^2) /
             cst.density_ice
        @test imp_solver.dtau_x[3, 2, 3] ≈ 2 * 0.9 / Λh
        # The Δτ field is now *uniform down the column*, which is the structural statement:
        # every trace of the sigma axis has left the explicit bound.
        @test all(≈(imp_solver.dtau_x[3, 2, 3]), interior(imp_solver.dtau_x))
    end

    @testset "one implicit step matches a reference tridiagonal solve" begin
        # The Thomas sweep is the one genuinely new piece of arithmetic in this phase, so it is
        # checked against `LinearAlgebra`'s own solve of the tridiagonal the *documentation*
        # specifies (source note above `_line_relax_x!`), assembled here from `µ`, `β`, `H` and
        # the sigma spacings rather than from the implementation. Uniform µ/H and dx = dy keep
        # the reference's `Λ_horiz` a closed form — every other coefficient is read off the
        # grid, including the non-uniform `QuadraticSigmaTransform` spacing.
        nx, ny, nz = 5, 4, 6
        dx, H0, μ0, β0 = 1e3, 500.0, 1e8, 2e3
        grid = StaggeredGrid(Float64, nx * dx, ny * dx, dx, dx, layering(; nz))
        rt = Runtime(grid)
        mech = MechanicState(grid)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0)
        fill_analytic!(mech.friction.beta_eff, rt.grid2d, (x, y) -> β0)

        momentum = BlatterPattynMomentumBalance()
        solver = PseudoTransientSolver(grid, momentum;
                                       vertical_treatment = ImplicitVertical(grid))
        pseudo_dt!(solver, mech, cst, rt, momentum)

        for k in 1:nz, j in -1:(ny + 2), i in -1:(nx + 2)
            solver.velocity_x_dt[i, j, k] = 0.1 * i + 0.01 * k - 0.03 * j
            solver.velocity_x_old[i, j, k] = 1.0 + 0.5 * k
            solver.velocity_y_dt[i, j, k] = -0.05 * j + 0.02 * k
            solver.velocity_y_old[i, j, k] = 2.0 - 0.25 * k
            mech.velocity.x[i, j, k] = solver.velocity_x_old[i, j, k]
            mech.velocity.y[i, j, k] = solver.velocity_y_old[i, j, k]
        end

        θ = 0.6
        Pagos._velocity_update!(solver.vertical_treatment, solver, mech, cst, rt, NoMask(),
                                mech.velocity.x, mech.velocity.y,
                                solver.velocity_x_old, solver.velocity_y_old, θ)

        Λh = 32μ0 / dx^2                                       # pre-mass-scaling, per §2.1
        τ̂  = 1 / Λh
        δζ(k) = Δz(rt.grid, Vertex(), 1, 1, k)
        Δζ(k) = Δz(rt.grid, Center(), 1, 1, k)
        A(k) = k > 1  ? μ0 / δζ(k) / Δζ(k) / H0^2 : 0.0        # bed-ward interface
        C(k) = k < nz ? μ0 / δζ(k + 1) / Δζ(k) / H0^2 : 0.0    # surface-ward interface
        D(k) = k == 1 ? β0 / Δζ(1) / H0 : 0.0                  # bed drag, diagonal only

        M = Tridiagonal([-τ̂ * A(k) for k in 2:nz],
                        [1 + τ̂ * (A(k) + C(k) + D(k)) for k in 1:nz],
                        [-τ̂ * C(k) for k in 1:(nz - 1)])
        for j in 1:ny, i in 1:(nx + 1)
            rhs = [θ * solver.dtau_x[i, j, k] * solver.velocity_x_dt[i, j, k] for k in 1:nz]
            expected = [solver.velocity_x_old[i, j, k] for k in 1:nz] .+ (M \ rhs)
            @test [mech.velocity.x[i, j, k] for k in 1:nz] ≈ expected rtol = 1e-12
        end
    end

    @testset "same fixed point as ExplicitVertical (masked, sheared, laterally varying)" begin
        nz = 6
        exp_case = sheared_margin_case(; nz)
        imp_case = sheared_margin_case(; nz)

        exp_res = solve_case(exp_case, ExplicitVertical())
        imp_res = solve_case(imp_case, ImplicitVertical(imp_case.grid))
        @test exp_res.res.converged
        @test imp_res.res.converged

        ux_e = interior(exp_case.mech.velocity.x)
        ux_i = interior(imp_case.mech.velocity.x)
        uy_e = interior(exp_case.mech.velocity.y)
        uy_i = interior(imp_case.mech.velocity.y)
        scale = maximum(abs, ux_e)

        @test scale > 1.0                                   # the case actually flows
        @test maximum(abs, ux_i .- ux_e) < 1e-2 * scale     # same fixed point, solver tol
        @test maximum(abs, uy_i .- uy_e) < 1e-2 * scale
        # Genuine vertical shear, or the comparison would be vacuous for this phase: over the
        # high-`β` band the column must be shear- rather than sliding-dominated.
        shear = maximum((abs(ux_e[i, j, nz] - ux_e[i, j, 1]) / max(abs(ux_e[i, j, 1]), 1e-6)
                         for i in axes(ux_e, 1), j in axes(ux_e, 2)))
        @test shear > 0.5
        # The margin column the strict-mask rule used to delete: finite and solved on both.
        @test all(isfinite, ux_i)
        @test isfinite(ux_i[exp_case.m + 1, 2, nz])
    end

    @testset "iteration count is independent of nz" begin
        # The whole claim of the phase. Explicit grows with the aspect ratio (`Λ_vert ∝ 1/Δz²`,
        # square-rooted by the damping); implicit is bounded by the horizontal operator alone
        # and so must be flat. Measured on the uniform slab rather than the masked case above
        # because a slab is the one geometry whose *solution* is nz-independent too, so any
        # growth is the solver's and not the discretization's.
        cc = (H0 = 1000.0, μ0 = 1e8, β0 = 1e3, α = 1e-2)
        dx, nx = 1e3, 8

        function slab_iterations(nz, vertical_treatment)
            grid = StaggeredGrid(Float64, nx * dx, 3dx, dx, dx, layering(; nz))
            rt, mech = Runtime(grid), MechanicState(grid)
            fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> cc.H0)
            fill_analytic!(mech.topography.surface, rt.grid2d, (x, y) -> -cc.α * x)
            fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> cc.μ0)
            fill_analytic!(mech.friction.beta_eff, rt.grid2d, (x, y) -> cc.β0)
            momentum = BlatterPattynMomentumBalance()
            solver = PseudoTransientSolver(grid, momentum;
                maxiter = 60_000, abstol = 1e-10, ncheck = 20,
                tuning = AutotunedDynamicRelaxation(),
                vertical_treatment = vertical_treatment(grid))
            res = pseudo_transient!(mech, cst, solver, rt, momentum)
            @test res.converged
            return res.iterations
        end

        nzs = (4, 8, 16, 32)
        expl = [slab_iterations(nz, _ -> ExplicitVertical()) for nz in nzs]
        impl = [slab_iterations(nz, ImplicitVertical) for nz in nzs]

        # Explicit pays the aspect ratio: an 8× refinement must cost several times more
        # iterations. Implicit must not — a factor of 2 of slack over the whole sweep, which
        # is far below the >10× explicit growth this is discriminating against.
        @test expl[end] > 4 * expl[1]
        @test maximum(impl) < 2 * minimum(impl)
        @test impl[end] < expl[end] / 4
    end
end

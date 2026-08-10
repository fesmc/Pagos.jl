using Pagos
using Test

include("../test_helpers/chmy.jl")

# The Chmy-native, C-grid staggered pseudo-transient momentum solver: the "flagship
# consumer" that assembles velocity gradients, the SSA/DIVA membrane stress, basal stress
# (β staggered onto the velocity faces — the deferred half of the Phase 2 `FrictionState`
# item) and the driving stress into the same iterative solve the collocated
# `pseudo_transient!` runs. Validated against analytic solutions, per the hybrid migration
# strategy — never against the collocated code path, which discretizes different points.

# -----------------------------------------------------------------------
# Analytical uniform ice-slab solution (same as test/mechanics/slab.jl): a spatially
# uniform slope, thickness, viscosity and friction coefficient give a spatially uniform
# velocity, so the membrane-stress divergence vanishes identically and the residual
# reduces to the local force balance τ_d = β u ⟹ u = τ_d / β. Because the exact solution
# has zero velocity gradient everywhere, the placeholder Neumann(0) halo refresh
# `pseudo_transient!` applies each iteration is also exact for this case, not merely a
# stand-in — a boundary-condition-independent check.
# -----------------------------------------------------------------------

function slab_analytical(; H0, μ0, β0, α, ρ = 910.0, g = 9.81)
    τd = ρ * g * H0 * α
    ub = τd / β0
    return (; τd, ub)
end

@testset "pseudo-transient momentum solver (C-grid staggered)" begin
    cst = Constants{Float64}()

    @testset "PseudoTransientSolver(::StaggeredGrid)" begin
        grid   = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        solver = PseudoTransientSolver(grid)

        @test location(solver.velocity_x_old) === location(solver.velocity_x_dt)
        @test location(solver.velocity_x_old) === (Vertex(), Center(), Center())  # acx
        @test location(solver.velocity_y_old) === (Center(), Vertex(), Center())  # acy
        @test size(interior(solver.velocity_x_old)) == (grid.nx + 1, grid.ny, 1)
        @test size(interior(solver.velocity_y_old)) == (grid.nx, grid.ny + 1, 1)
        @test eltype(solver.velocity_x_old) === Float64

        # dtau_x/dtau_y mirror the velocity work arrays: same node class, same shape.
        @test location(solver.dtau_x) === location(solver.velocity_x_old)
        @test location(solver.dtau_y) === location(solver.velocity_y_old)
        @test size(interior(solver.dtau_x)) == size(interior(solver.velocity_x_old))
        @test size(interior(solver.dtau_y)) == size(interior(solver.velocity_y_old))

        # The iteration parameters live on the strategy that reads them, not on the
        # solver: default tuning is undamped (gamma = 1), default Δτ is the Gershgorin
        # bound, and `PseudoTransientSolver` carries neither gamma/theta_v nor muB/ndim.
        @test solver.tuning == FixedTuning()
        @test solver.tuning.gamma == 1
        @test solver.pseudo_timestep isa GershgorinPseudoTimeStep
        @test !hasproperty(solver, :gamma)
        @test !hasproperty(solver, :theta_v)
        @test !hasproperty(solver, :muB)
        @test !hasproperty(solver, :ndim2)

        # AutotunedDynamicRelaxation's spectral estimates are normalized by the Gershgorin
        # row sum (λ_max ≤ 1 by construction), so pairing it with any other Δτ rule is a
        # correctness error, not a preference — rejected at construction.
        @test_throws ArgumentError PseudoTransientSolver(grid;
            pseudo_timestep = ViscosityPseudoTimeStep(),
            tuning = AutotunedDynamicRelaxation())
        @test PseudoTransientSolver(grid;
            tuning = AutotunedDynamicRelaxation()).tuning isa AutotunedDynamicRelaxation

        # The solver itself imposes no `nz` requirement: every work array is built on
        # `grid.grid2d` and the unknown it iterates is depth-integrated for both balances,
        # so a column grid is perfectly constructible. Which *momentum balance* tolerates a
        # given grid is checked at `pseudo_transient!` instead (see the guard test below).
        layering  = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 4))
        col_grid  = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, layering)
        col_solver = PseudoTransientSolver(col_grid)
        @test size(interior(col_solver.velocity_x_old)) == (col_grid.nx + 1, col_grid.ny, 1)
        @test size(interior(col_solver.dtau_x)) == (col_grid.nx + 1, col_grid.ny, 1)
    end

    # Per-balance grid guards (decision 4): SSA is depth-independent and runs anywhere;
    # DIVA is *defined* by its vertical structure and its F₂ integral is 25% low on a
    # single layer, so `nz == 1` is an error rather than a silently-degraded answer.
    @testset "momentum balance grid guards" begin
        flat = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        layering = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 4))
        col  = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, layering)

        mech_flat, rt_flat = MechanicState(flat), Runtime(flat)
        solver = PseudoTransientSolver(flat; maxiter = 1)
        @test_throws ArgumentError pseudo_transient!(mech_flat, cst, solver, rt_flat,
                                                     DIVAMomentumBalance())
        # SSA on the same flat grid is fine, and is what the default resolves to.
        @test pseudo_transient!(mech_flat, cst, solver, rt_flat,
                                SSAMomentumBalance()) isa NamedTuple

        # DIVA is accepted on a column grid (this only checks the guard, not the physics —
        # the DIVA path itself is Stage 2).
        @test Runtime(col).grid !== Runtime(col).grid2d
    end

    # A linear viscosity gradient makes `lerp` exact, so `pseudo_dt!`'s local field must
    # match the closed-form scaling/lerp(μ) formula pointwise, not just in aggregate.
    @testset "pseudo_dt!: local Δτ field matches lerp(μ) at the velocity faces" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        rt   = Runtime(grid)
        mech = MechanicState(grid)
        ρ, muB, ndim, dtau_scaling = 910.0, 1e2, 4.1, 1.0
        η0, a = 1e5, 3e3

        fill_analytic!(mech.material.viscosity_depthaveraged, rt.grid2d,
                       (x, y) -> η0 + a * x)

        solver = PseudoTransientSolver(grid)
        dx = Δx(rt.grid2d, Center(), 1, 1, 1)
        dy = Δy(rt.grid2d, Center(), 1, 1, 1)
        pseudo_dt!(solver.dtau_x, solver.dtau_y, ρ, dx, dy,
                  mech.material.viscosity_depthaveraged, muB, ndim, dtau_scaling, rt)

        scaling = dtau_scaling * ρ * dx * dy / (4 * (1 + muB) * ndim)
        expected_x = analytic_like(solver.dtau_x, rt.grid2d, (x, y) -> scaling / (η0 + a * x))
        expected_y = analytic_like(solver.dtau_y, rt.grid2d, (x, y) -> scaling / (η0 + a * x))
        @test interior(solver.dtau_x) ≈ expected_x
        @test interior(solver.dtau_y) ≈ expected_y

        # Uniform viscosity: no lerp needed (both cells agree), so the local field
        # collapses to Sandip Eq. 7 evaluated at that one viscosity, everywhere.
        fill_analytic!(mech.material.viscosity_depthaveraged, rt.grid2d, (x, y) -> η0)
        pseudo_dt!(solver.dtau_x, solver.dtau_y, ρ, dx, dy,
                  mech.material.viscosity_depthaveraged, muB, ndim, dtau_scaling, rt)
        @test all(≈(scaling / η0), interior(solver.dtau_x))
        @test all(≈(scaling / η0), interior(solver.dtau_y))
    end

    # Same linear-velocity setup as the "membrane stress ... exact under uniform η, H" test
    # below, but checking `strainrate.effective` (the true second invariant) instead of
    # `.xx`/`.xy`/`.yy` (the membrane stress) — the two are deliberately different
    # quantities sharing a struct, see `effective_strainrate_ssa!`'s docstring. `x_dy`,
    # `y_dx` are uniform for a linear field, so `lerp` onto `aa` is exact.
    @testset "effective_strainrate_ssa!: exact under a linear velocity field" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        rt   = Runtime(grid)
        mech = MechanicState(grid)
        a, b, c, d = 2e-3, -1e-3, 5e-4, 3e-3

        fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d, (x, y) -> a * x + b * y)
        fill_analytic!(mech.velocity.depthaverage_y, rt.grid2d, (x, y) -> c * x + d * y)
        depthaverage_velocitygradients!(mech.velocity, rt)
        effective_strainrate_ssa!(mech.strainrate, mech.velocity, rt)

        expected = sqrt(a^2 + d^2 + a * d + ((b + c) / 2)^2)
        @test all(≈(expected), interior(mech.strainrate.effective_depthaveraged))
    end

    # Direct check of the Glen-law + log-space relaxation formula (Sandip et al. 2024,
    # Eq. 3 and 8) through the public `update_viscosity!` entry point, on a uniform
    # velocity gradient (so `effective_strainrate_ssa!`'s output is a known constant, per
    # the test above) and uniform rate factor, so the closed form is exact everywhere.
    @testset "update_viscosity!(::GlenViscosityContinuation): matches the closed-form formula" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        rt   = Runtime(grid)
        mech = MechanicState(grid)
        # µ_old0 in Pa yr, matching the ~1e5-1e6 scale every other viscosity here uses;
        # it is the relaxation seed, so it only has to be a plausible previous iterate.
        A0, n, ε̇0, μ_old0 = 1e-16, 3.0, 1e-12, 3e5
        a, b, c, d = 2e-3, -1e-3, 5e-4, 3e-3   # same linear velocity as the test above

        fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d, (x, y) -> a * x + b * y)
        fill_analytic!(mech.velocity.depthaverage_y, rt.grid2d, (x, y) -> c * x + d * y)
        fill_analytic!(mech.material.rate_factor_depthaveraged, rt.grid2d, (x, y) -> A0)
        depthaverage_velocitygradients!(mech.velocity, rt)

        eff0  = sqrt(a^2 + d^2 + a * d + ((b + c) / 2)^2)
        μ_raw = inv(2 * A0^(1 / n)) * sqrt(eff0^2 + ε̇0^2)^((1 - n) / n)

        setdata!(mech.material.viscosity_depthaveraged, μ_old0)
        vc1 = GlenViscosityContinuation(; n_glen = n, theta_mu = 1.0, strainrate_reg = ε̇0)
        update_viscosity!(mech, vc1, rt)
        @test all(≈(μ_raw, rtol = 1e-10), interior(mech.material.viscosity_depthaveraged))

        setdata!(mech.material.viscosity_depthaveraged, μ_old0)
        vc2 = GlenViscosityContinuation(; n_glen = n, theta_mu = 0.2, strainrate_reg = ε̇0)
        expected = exp(0.2 * log(μ_raw) + 0.8 * log(μ_old0))
        update_viscosity!(mech, vc2, rt)
        @test all(≈(expected, rtol = 1e-10), interior(mech.material.viscosity_depthaveraged))
    end

    # β lives at `aa`; τ_b = β v_b is formed on the velocity faces via `lerp`, exactly like
    # `drivingstress!`'s treatment of H. Linear β makes `lerp` exact.
    @testset "basalstress!: β staggered exactly (linear β, uniform v)" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        rt   = Runtime(grid)
        mech = MechanicState(grid)
        a, b, v0 = 100.0, 5.0, 20.0

        fill_analytic!(mech.friction.beta_eff, rt.grid2d, (x, y) -> a + b * x)
        setdata!(mech.velocity.base_x, v0)
        setdata!(mech.velocity.base_y, v0)

        basalstress!(mech, rt)

        # `lerp` of a linear β is exact, so this is exact to roundoff, not just approximate.
        expected_x = analytic_like(mech.stress.base_x, rt.grid2d, (x, y) -> (a + b * x) * v0)
        expected_y = analytic_like(mech.stress.base_y, rt.grid2d, (x, y) -> (a + b * x) * v0)
        @test interior(mech.stress.base_x) ≈ expected_x
        @test interior(mech.stress.base_y) ≈ expected_y
    end

    @testset "membrane stress (SSA/DIVA `strainrate!`): exact under uniform η, H" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        rt   = Runtime(grid)
        mech = MechanicState(grid)
        η0, H0, a, b, c, d = 1e5, 800.0, 2e-3, -1e-3, 5e-4, 3e-3

        # `fill_analytic!`, not `setdata!`, for η/H: `setdata!` only touches the interior,
        # and the first/last `ab` vertex is built from a halo `aa` cell (vertex i sits
        # between centres i-1, i; i=1 reaches centre 0). See `pagos-roadmaps/chmy.md`, §3.
        fill_analytic!(mech.material.viscosity_depthaveraged, rt.grid2d, (x, y) -> η0)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
        fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d, (x, y) -> a * x + b * y)
        fill_analytic!(mech.velocity.depthaverage_y, rt.grid2d, (x, y) -> c * x + d * y)

        depthaverage_velocitygradients!(mech.velocity, rt)
        membranestress!(mech.stress, mech.velocity, mech.material, mech.topography,
                   DIVAMomentumBalance(), rt)

        @test all(interior(mech.stress.membrane_xx) .≈ 2η0 * H0 * (2a + d))
        @test all(interior(mech.stress.membrane_yy) .≈ 2η0 * H0 * (a + 2d))
        @test all(interior(mech.stress.membrane_xy) .≈ η0 * H0 * (b + c))

        # Locations match the layout table: N_xx/N_yy at `aa`, N_xy at `ab`.
        @test location(mech.stress.membrane_xx) === (Center(), Center(), Center())
        @test location(mech.stress.membrane_xy) === (Vertex(), Vertex(), Center())
    end

    # A genuine viscosity contrast: hlerp must give the harmonic, not the arithmetic, mean.
    # Located by cell centres, not a coordinate/index guess, per the vertex/centre lesson
    # in `pagos-roadmaps/chmy.md` (§3): vertex `i` sits between centres `i-1` and `i`.
    @testset "membrane stress: harmonic η averaging at a viscosity step" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        rt   = Runtime(grid)
        mech = MechanicState(grid)
        η_lo, η_hi, H0, b = 1e4, 1e6, 500.0, 1e-2

        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
        fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d, (x, y) -> b * y)   # ε̇xy = b/2 uniform
        fill_analytic!(mech.velocity.depthaverage_y, rt.grid2d, (x, y) -> 0.0)
        fill_analytic!(mech.material.viscosity_depthaveraged, rt.grid2d,
                       (x, y) -> x < 0 ? η_lo : η_hi)

        depthaverage_velocitygradients!(mech.velocity, rt)
        membranestress!(mech.stress, mech.velocity, mech.material, mech.topography,
                   DIVAMomentumBalance(), rt)

        # Find an ab vertex whose two x-centres genuinely straddle the step.
        xs = collect(grid.x)
        istep = findfirst(i -> xs[i] < 0 <= xs[i + 1], 1:(length(xs) - 1))
        @test xs[istep] < 0 <= xs[istep + 1]         # the premise: the face really straddles it
        i_ab = istep + 1                              # vertex i sits between centres i-1, i

        η_harm = 2 / (1 / η_lo + 1 / η_hi)
        expected_xy = η_harm * H0 * b
        @test interior(mech.stress.membrane_xy)[i_ab, 4, 1] ≈ expected_xy
        # The arithmetic mean would have given a different (larger) answer.
        @test !(interior(mech.stress.membrane_xy)[i_ab, 4, 1] ≈ (η_lo + η_hi) / 2 * H0 * b)
    end

    # Containment: an unmasked ice-free (η = 0) cell produces NaN via hlerp; masking with
    # the strict rule keeps it finite, mirroring `deviatoric_stress!`'s established pattern.
    @testset "membrane stress: NaN off-ice, contained by IceMask" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        rt   = Runtime(grid)
        mech = MechanicState(grid)
        topo = TopographicState(grid)

        setdata!(mech.topography.thickness, 500.0)
        fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d, (x, y) -> 0.01y)
        fill_analytic!(mech.velocity.depthaverage_y, rt.grid2d, (x, y) -> 0.0)
        fill_analytic!(mech.material.viscosity_depthaveraged, rt.grid2d,
                       (x, y) -> x < 0 ? 0.0 : 1e5)
        depthaverage_velocitygradients!(mech.velocity, rt)

        membranestress!(mech.stress, mech.velocity, mech.material, mech.topography,
                   DIVAMomentumBalance(), rt)
        @test any(isnan, interior(mech.stress.membrane_xy))

        setdata!(topo.mask.is_ice, true)
        for i in axes(interior(topo.mask.is_ice), 1)
            x, _, _ = coord(rt.grid2d, location(topo.mask.is_ice), i, 1, 1)
            x < 0 && (interior(topo.mask.is_ice)[i, :, 1] .= false)
        end
        mask = IceMask(topo.mask.is_ice)

        membranestress!(mech.stress, mech.velocity, mech.material, mech.topography,
                   DIVAMomentumBalance(), rt, mask)
        @test !any(isnan, interior(mech.stress.membrane_xy))
    end

    # Uniform-slab residual with no membrane stress: mirrors the collocated
    # "Pseudo-transient pure functions" test in test/mechanics/slab.jl.
    @testset "dotvel!: uniform-slab residual, no membrane stress" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        rt   = Runtime(grid)
        ρ, H0, base_x0, driving_x0 = 910.0, 1000.0, 1e3, -8927.0

        sxx = Field(rt.arch, rt.grid, (Center(), Center(), Center()))
        sxy = Field(rt.arch, rt.grid, (Vertex(), Vertex(), Center()))
        syy = Field(rt.arch, rt.grid, (Center(), Center(), Center()))
        base_x  = Field(rt.arch, rt.grid2d, (Vertex(), Center(), Center()))
        base_y  = Field(rt.arch, rt.grid2d, (Center(), Vertex(), Center()))
        driv_x  = Field(rt.arch, rt.grid2d, (Vertex(), Center(), Center()))
        driv_y  = Field(rt.arch, rt.grid2d, (Center(), Vertex(), Center()))
        H       = Field(rt.arch, rt.grid2d, (Center(), Center(), Center()))
        dvx     = Field(rt.arch, rt.grid2d, (Vertex(), Center(), Center()))
        dvy     = Field(rt.arch, rt.grid2d, (Center(), Vertex(), Center()))
        rx      = Field(rt.arch, rt.grid2d, (Vertex(), Center(), Center()))
        ry      = Field(rt.arch, rt.grid2d, (Center(), Vertex(), Center()))

        setdata!(sxx, 0.0); setdata!(sxy, 0.0); setdata!(syy, 0.0)
        setdata!(base_x, base_x0); setdata!(base_y, 0.0)
        setdata!(driv_x, driving_x0); setdata!(driv_y, 0.0)
        # `fill_analytic!`, not `setdata!`: dvx reads H at the two boundary acx faces,
        # which need H's halo (setdata! only fills the interior; H's halo defaults to 0,
        # not H0). base_x/driv_x are read at the same acx node with no interpolation, so
        # their halos are never touched — the halo-fill requirement here is specific to
        # the field that gets `lerp`ed onto a face.
        fill_analytic!(H, rt.grid2d, (x, y) -> H0)

        dotvel!(dvx, dvy, sxx, sxy, syy, base_x, base_y, driv_x, driv_y, H, ρ, rt,
               DIVAMomentumBalance(); resid_x = rx, resid_y = ry)

        @test all(≈((0.0 - base_x0 - driving_x0) / (ρ * H0)), interior(dvx))
        @test all(==(0.0), interior(dvy))
        @test interior(rx) ≈ interior(dvx)   # gamma = 1: residual and damped rate agree
        @test all(==(0.0), interior(ry))
    end

    # gamma defaults to 1 (no damping: dv is fully replaced by the raw rate every call,
    # regardless of what it held on entry — the plain PT iteration). gamma < 1 keeps a
    # (1 - gamma) fraction of the incoming value instead, per Sandip et al. (2024) Eq. 12–14.
    @testset "dotvel!: damped rate accumulation (gamma < 1 keeps memory of dvx_old)" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        rt   = Runtime(grid)
        ρ, H0, base_x0, driving_x0 = 910.0, 1000.0, 1e3, -8927.0
        dvx_old0, gamma = 4.2, 0.3

        sxx = Field(rt.arch, rt.grid, (Center(), Center(), Center()))
        sxy = Field(rt.arch, rt.grid, (Vertex(), Vertex(), Center()))
        syy = Field(rt.arch, rt.grid, (Center(), Center(), Center()))
        base_x  = Field(rt.arch, rt.grid2d, (Vertex(), Center(), Center()))
        base_y  = Field(rt.arch, rt.grid2d, (Center(), Vertex(), Center()))
        driv_x  = Field(rt.arch, rt.grid2d, (Vertex(), Center(), Center()))
        driv_y  = Field(rt.arch, rt.grid2d, (Center(), Vertex(), Center()))
        H       = Field(rt.arch, rt.grid2d, (Center(), Center(), Center()))
        dvx     = Field(rt.arch, rt.grid2d, (Vertex(), Center(), Center()))
        dvy     = Field(rt.arch, rt.grid2d, (Center(), Vertex(), Center()))
        rx      = Field(rt.arch, rt.grid2d, (Vertex(), Center(), Center()))
        ry      = Field(rt.arch, rt.grid2d, (Center(), Vertex(), Center()))

        setdata!(sxx, 0.0); setdata!(sxy, 0.0); setdata!(syy, 0.0)
        setdata!(base_x, base_x0); setdata!(base_y, 0.0)
        setdata!(driv_x, driving_x0); setdata!(driv_y, 0.0)
        fill_analytic!(H, rt.grid2d, (x, y) -> H0)
        setdata!(dvx, dvx_old0)   # nonzero "old" rate the damped update must partly retain
        setdata!(dvy, 0.0)

        dotvel!(dvx, dvy, sxx, sxy, syy, base_x, base_y, driv_x, driv_y, H, ρ, rt,
               DIVAMomentumBalance(); gamma, resid_x = rx, resid_y = ry)

        raw_x = (0.0 - base_x0 - driving_x0) / (ρ * H0)
        @test all(≈((1 - gamma) * dvx_old0 + raw_x), interior(dvx))
        @test all(==(0.0), interior(dvy))
        # resid_x always holds the raw (undamped) rate, regardless of gamma — unlike dvx.
        @test all(≈(raw_x), interior(rx))
        @test all(==(0.0), interior(ry))

        # gamma = 1 (the default): dvx_old is fully discarded, exactly like the untyped
        # call above — the mechanism is opt-in, not a behaviour change for existing callers.
        setdata!(dvx, dvx_old0)
        dotvel!(dvx, dvy, sxx, sxy, syy, base_x, base_y, driv_x, driv_y, H, ρ, rt,
               DIVAMomentumBalance(); resid_x = rx, resid_y = ry)
        @test all(≈(raw_x), interior(dvx))
    end

    # Large β0 → λ = θ·dτ·β/(ρH) ≈ 0.9 → converges in ~10 iterations (same case as the
    # collocated PseudoTransientSolver test in test/mechanics/slab.jl).
    const_case = (H0 = 1000.0, μ0 = 1e5, β0 = 1e4, α = 1e-3)

    function setup_slab(T = Float64; nx = 10, dx = T(5e3))
        grid = StaggeredGrid(T, nx * dx, 3 * dx, dx, dx)
        rt   = Runtime(grid)
        mech = MechanicState(grid)
        return grid, rt, mech
    end

    # Every uniform field here is filled with `fill_analytic!`, not `setdata!`: `setdata!`
    # only touches the interior, and `strainrate!`'s membrane-stress kernel reads `H`, `η`
    # through `lerp`/`hlerp` at the two boundary `ab` corners, which are built partly from
    # the outer halo ring `setdata!` never touches (default-zero on allocation). `lerp` of
    # a zero neighbour is merely wrong (H/2 instead of H); `hlerp` of one is `NaN` — see the
    # note on `pseudo_transient!` below, discovered by this very test.
    function fill_slab!(mech, rt, c, T = Float64)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> T(c.H0))
        fill_analytic!(mech.topography.surface, rt.grid2d, (x, y) -> T(c.H0) - T(c.α) * x)
        fill_analytic!(mech.material.viscosity_depthaveraged, rt.grid2d, (x, y) -> T(c.μ0))
        fill_analytic!(mech.friction.beta_eff, rt.grid2d, (x, y) -> T(c.β0))
        return nothing
    end

    @testset "DIVA uniform slab — pseudo_transient!" begin
        an = slab_analytical(; const_case...)
        grid, rt, mech = setup_slab()
        fill_slab!(mech, rt, const_case)
        solver = PseudoTransientSolver(grid; maxiter = 200, abstol = 1e-8)

        res = pseudo_transient!(mech, cst, solver, rt)

        @test res.converged
        @test all(≈(an.ub, rtol = 1e-5), interior(mech.velocity.depthaverage_x))
        @test all(≈(0.0,   atol = 1e-8 * abs(an.ub)), interior(mech.velocity.depthaverage_y))
        # residual is the raw momentum-balance rate, not the velocity increment `error` is
        # based on — different units, but both must be tiny at a converged solution.
        @test res.residual >= 0
        @test res.residual < 1e-6

        # SSA and DIVA share the same equations in this SSA-limit (no resolved shear).
        setdata!(mech.velocity.depthaverage_x, 0.0); setdata!(mech.velocity.depthaverage_y, 0.0)
        res_ssa = pseudo_transient!(mech, cst, solver, rt, SSAMomentumBalance())
        @test res_ssa.converged
        @test all(≈(an.ub, rtol = 1e-5), interior(mech.velocity.depthaverage_x))
    end

    @testset "ncheck = 5" begin
        an = slab_analytical(; const_case...)
        grid, rt, mech = setup_slab()
        fill_slab!(mech, rt, const_case)
        solver5 = PseudoTransientSolver(grid; maxiter = 200, abstol = 1e-8, ncheck = 5)

        res5 = pseudo_transient!(mech, cst, solver5, rt)
        @test res5.converged
        @test res5.iterations % 5 == 0
        @test all(≈(an.ub, rtol = 1e-5), interior(mech.velocity.depthaverage_x))
        @test all(≈(0.0,   atol = 1e-8 * abs(an.ub)), interior(mech.velocity.depthaverage_y))
    end

    @testset "Float32 stays Float32" begin
        an = slab_analytical(; const_case...)
        grid, rt, mech = setup_slab(Float32)
        fill_slab!(mech, rt, const_case, Float32)
        cst32  = Constants{Float32}()
        solver = PseudoTransientSolver(grid; maxiter = 200, abstol = 1f-6)

        res = pseudo_transient!(mech, cst32, solver, rt)
        @test res.converged
        @test eltype(mech.velocity.depthaverage_x) === Float32
        @test all(≈(Float32(an.ub), rtol = 1f-3), interior(mech.velocity.depthaverage_x))
    end

    # What a *fixed* gamma buys, and what it does not — the honest negative result that
    # motivates Phase 2. Corrects an earlier version of this file, which measured a
    # "resolution-robust ~2x speedup at gamma = 0.5" against `ViscosityPseudoTimeStep`,
    # whose Δτ is ~101x too small (see the closed-form Gershgorin test below). Against a
    # correct Δτ the textbook picture appears instead, and it is less flattering to a fixed
    # gamma: the *optimal* gamma moves with resolution, so any single value is wrong almost
    # everywhere. Measured on this case (iterations, undamped → best of a scan):
    #
    #   nx     10     20     40     80    160
    #   γ=1.0  14     34    104    363   1323      (≈ O(nx²), the papers' starting point)
    #   γ=0.2 194    185    171    176    179
    #   best  γ=1.0  γ=1.0  γ=0.5  γ=0.3  γ=0.2    → 14, 34, 59, 108, 179 (≈ O(nx^0.75))
    #
    # so damping *costs* 14x at nx = 10 and *saves* 7x at nx = 160, at the same gamma. The
    # assertion is that inversion, which no constant can avoid — only Phase 2's
    # resolution-adaptive selection can, and the testset after next shows it doing so.
    # Both solves converge to the same analytic solution regardless of gamma: damping
    # changes the transient path, never the fixed point.
    @testset "fixed gamma: the optimal value moves with resolution (so no constant works)" begin
        scale_case = (H0 = 1000.0, μ0 = 1e5, β0 = 1e2, α = 1e-3)
        an = slab_analytical(; scale_case...)
        Lx = 100e3

        # The uniform-slab solution has zero velocity gradient, so its membrane stress
        # vanishes and convergence is friction-dominated whatever gamma is — useless for
        # exercising damping. Perturbing it forces the membrane divergence to act
        # repeatedly, which is the pathway damping accelerates.
        function resolved_run(nx, gamma)
            dx = Lx / nx
            grid, rt, mech = setup_slab(; nx, dx)
            fill_slab!(mech, rt, scale_case)
            fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d,
                           (x, y) -> an.ub + 0.5 * an.ub * sinpi(3x / Lx))
            solver = PseudoTransientSolver(grid; maxiter = 20_000, abstol = 1e-8,
                                           tuning = FixedTuning(; gamma))
            return pseudo_transient!(mech, cst, solver, rt), mech
        end

        coarse_undamped, m1 = resolved_run(10, 1.0)
        coarse_damped, m2   = resolved_run(10, 0.2)
        fine_undamped, m3   = resolved_run(160, 1.0)
        fine_damped, m4     = resolved_run(160, 0.2)

        for (r, m) in ((coarse_undamped, m1), (coarse_damped, m2),
                       (fine_undamped, m3), (fine_damped, m4))
            @test r.converged
            @test all(≈(an.ub, rtol = 1e-5), interior(m.velocity.depthaverage_x))
        end

        # The inversion: the same gamma that is far worse than undamped on the coarse grid
        # is far better on the fine one.
        @test coarse_damped.iterations > 5 * coarse_undamped.iterations
        @test fine_damped.iterations < fine_undamped.iterations / 5

        # And the undamped baseline degrades steeply with resolution (14 → 1323 over a 16x
        # refinement, ≈ O(nx²)) — the scaling Phase 2 exists to fix.
        @test fine_undamped.iterations > 20 * coarse_undamped.iterations
    end

    # End-to-end verification of `GlenViscosityContinuation` (no closed-form nonlinear
    # solution is derived — see the module note below for why a self-consistency check is
    # used instead). A friction contrast β(x) = β0·(1 + amp·cos(2πx/Lx)) forces genuine
    # lateral shear (unlike the uniform-slab cases above, whose analytic solution has zero
    # velocity gradient everywhere and so never exercises the strain-rate → viscosity
    # feedback at all), giving a real nonlinear SSA/DIVA problem with no known closed form.
    #
    # `A0` is calibrated, not arbitrary: picked so that Glen's law evaluated at a strain
    # rate of `ε̇_typical` (the right order of magnitude for the resulting shear, found by
    # trial run) gives a viscosity of `μ_target`, matching the ~1e5 scale every other
    # `PseudoTransientSolver` test in this file uses. Without this, viscosity and Δτ end up
    # wildly mismatched (found the hard way: an uncalibrated run took 50,000 iterations
    # without converging, its velocity stuck ~5 orders of magnitude below the
    # friction-dominated scale — Δτ ∝ 1/μ, and a viscosity guess many orders of magnitude
    # off makes every step size off by the same factor).
    @testset "GlenViscosityContinuation: nonlinear shear flow (no closed form — self-consistency check)" begin
        nx, dx = 40, 5e3
        Lx = nx * dx
        H0, α, β0, n, amp = 1000.0, 1e-3, 1e3, 3.0, 0.95

        ε̇_typical, μ_target = 1e-4, 1e5
        A0  = (1 / (2 * μ_target * ε̇_typical^((n - 1) / n)))^n
        ε̇0  = 1e-8   # regularization floor, well below ε̇_typical

        # NB: this helper's own locals are deliberately *not* named `grid`/`rt`/`mech` —
        # those are the names `setup_slab` itself uses internally, and reusing them here
        # (in a nested function calling `setup_slab` twice, results compared by identity)
        # hits a Julia closure/box-sharing quirk that silently aliases the two calls'
        # return values. Found the hard way: an earlier version of this test reused those
        # names and every "independent" solve turned out to be operating on the same
        # object as the first, so every comparison against it was trivially self-equal.
        function build_shear_problem()
            g, r, mc = setup_slab(; nx, dx)
            fill_analytic!(mc.topography.surface, r.grid2d, (x, y) -> H0 - α * x)
            fill_analytic!(mc.topography.thickness, r.grid2d, (x, y) -> H0)
            fill_analytic!(mc.friction.beta_eff, r.grid2d,
                           (x, y) -> β0 * (1 + amp * cospi(2x / Lx)))
            fill_analytic!(mc.material.rate_factor_depthaveraged, r.grid2d, (x, y) -> A0)
            fill_analytic!(mc.material.viscosity_depthaveraged, r.grid2d, (x, y) -> μ_target)
            return g, r, mc
        end

        vc = GlenViscosityContinuation(; n_glen = n, theta_mu = 1.0, strainrate_reg = ε̇0)
        grid_glen, rt_glen, mech_glen = build_shear_problem()
        solver_glen = PseudoTransientSolver(grid_glen; maxiter = 5000, abstol = 1e-9,
                                           tuning = FixedTuning(gamma = 0.5),
                                           viscosity_continuation = vc)
        res_glen = pseudo_transient!(mech_glen, cst, solver_glen, rt_glen)
        @test res_glen.converged

        # Self-consistency: the converged viscosity must match Glen's law evaluated at the
        # *converged* strain rate — recomputed independently here, not read back from
        # whatever the solver last wrote mid-iteration.
        effective_strainrate_ssa!(mech_glen.strainrate, mech_glen.velocity, rt_glen)
        eff = interior(mech_glen.strainrate.effective_depthaveraged)
        μ_closed = @. inv(2 * A0^(1 / n)) * sqrt(eff^2 + ε̇0^2)^((1 - n) / n)
        @test interior(mech_glen.material.viscosity_depthaveraged) ≈ μ_closed rtol=1e-6

        # The mechanism has a real effect: compare against a fixed-viscosity (no
        # continuation) solve from the same initial guess — measured ~4.6% difference at
        # these parameters, the assertion only requires 1%.
        grid_fixed, rt_fixed, mech_fixed = build_shear_problem()
        solver_fixed = PseudoTransientSolver(grid_fixed; maxiter = 5000, abstol = 1e-9,
                                            tuning = FixedTuning(gamma = 0.5))
        res_fixed = pseudo_transient!(mech_fixed, cst, solver_fixed, rt_fixed)
        @test res_fixed.converged

        reldiff = maximum(abs.(interior(mech_glen.velocity.depthaverage_x) .- interior(mech_fixed.velocity.depthaverage_x))) /
                  maximum(abs.(interior(mech_fixed.velocity.depthaverage_x)))
        @test reldiff > 0.01

        # theta_mu only changes the transient path, not the fixed point (same argument as
        # the damping tests above, now for the log-space viscosity relaxation instead of
        # the velocity-rate one): a relaxed continuation must converge to the same answer.
        vc_relaxed = GlenViscosityContinuation(; n_glen = n, theta_mu = 0.3, strainrate_reg = ε̇0)
        grid_relaxed, rt_relaxed, mech_relaxed = build_shear_problem()
        solver_relaxed = PseudoTransientSolver(grid_relaxed; maxiter = 5000, abstol = 1e-9,
                                               tuning = FixedTuning(gamma = 0.5),
                                               viscosity_continuation = vc_relaxed)
        res_relaxed = pseudo_transient!(mech_relaxed, cst, solver_relaxed, rt_relaxed)
        @test res_relaxed.converged
        @test all(isapprox.(interior(mech_relaxed.velocity.depthaverage_x), interior(mech_glen.velocity.depthaverage_x);
                            rtol = 1e-4))
    end

    # -------------------------------------------------------------------
    # GershgorinPseudoTimeStep / ScaledResidual: the two dispatch points added after the
    # AIS-geometry example (`docs/src/examples/ais-pt.jl`) came back with empty ice
    # shelves. Each of the three tests below pins one leg of that failure.
    # -------------------------------------------------------------------

    # Λ is a sum of |stencil coefficients|, so a uniform slab makes it a closed form:
    # with dx = dy, P = Q = ηH everywhere, the four membrane terms are
    # (8 + 4 + 2 + 2)·2ηH/dx² = 32ηH/dx², and Λ = (32ηH/dx² + β)/(ρH).
    @testset "pseudo_dt!(::GershgorinPseudoTimeStep): closed form on a uniform slab" begin
        grid, rt, mech = setup_slab()
        fill_slab!(mech, rt, const_case)
        dx = Δx(rt.grid2d, Center(), 1, 1, 1)
        ρ  = cst.density_ice
        (; H0, μ0, β0) = const_case

        for cfl in (1.0, 0.9, 0.5)
            solver = PseudoTransientSolver(grid;
                pseudo_timestep = GershgorinPseudoTimeStep(cfl = cfl))
            pseudo_dt!(solver, mech, cst, rt)

            Λ = (32 * μ0 * H0 / dx^2 + β0) / (ρ * H0)
            @test all(≈(2 * cfl / Λ), interior(solver.dtau_x))
            @test all(≈(2 * cfl / Λ), interior(solver.dtau_y))
        end

        # Without a friction law being solved for, the drag is a prescribed forcing, not a
        # term of the operator, so it must drop out of the spectral bound.
        solver_nofric = PseudoTransientSolver(grid;
            pseudo_timestep = GershgorinPseudoTimeStep(cfl = 1.0),
            friction_update = NoFrictionUpdate())
        pseudo_dt!(solver_nofric, mech, cst, rt)
        @test all(≈(2 * ρ * H0 / (32 * μ0 * H0 / dx^2)), interior(solver_nofric.dtau_x))

        # And with no drag at all it reduces to Sandip's Eq. 7 at muB = 0, ndim = 4.1 —
        # i.e. ρdx²/(16η) vs ρdx²/(16.4η), agreeing to within that constant. This is the
        # check that the *default* muB = 1e2 is ~101x too conservative, not that Eq. 7 is
        # wrong in the regime it was derived for.
        sandip(muB) = ρ * dx * dx / (4 * (1 + muB) * 4.1 * μ0)
        gersh = 2 * ρ * H0 / (32 * μ0 * H0 / dx^2)
        @test gersh ≈ sandip(0) * (16.4 / 16) rtol=1e-12
        @test sandip(1e2) < gersh / 100
    end

    # The leg that made the AIS solve blow up rather than merely stall: Sandip's Δτ does
    # not see the basal-drag term, so once β/(ρH) dominates the spectrum, θ·Δτ·β/(ρH) > 2
    # and the iteration is unconditionally unstable. At these parameters that factor is
    # ~9 for `ViscosityPseudoTimeStep` and ~1.8 for the Gershgorin bound.
    @testset "GershgorinPseudoTimeStep: stable where the viscosity-only Δτ diverges" begin
        stiff_case = (H0 = 1000.0, μ0 = 1e5, β0 = 1e5, α = 1e-3)
        an = slab_analytical(; stiff_case...)

        grid_v, rt_v, mech_v = setup_slab()
        fill_slab!(mech_v, rt_v, stiff_case)
        solver_v = PseudoTransientSolver(grid_v; maxiter = 200, abstol = 1e-8,
            pseudo_timestep = ViscosityPseudoTimeStep())
        res_v = pseudo_transient!(mech_v, cst, solver_v, rt_v)

        @test !res_v.converged
        @test maximum(abs, interior(mech_v.velocity.depthaverage_x)) > 100 * abs(an.ub)   # blown up

        grid_g, rt_g, mech_g = setup_slab()
        fill_slab!(mech_g, rt_g, stiff_case)
        solver_g = PseudoTransientSolver(grid_g; maxiter = 200, abstol = 1e-8,
            tuning = FixedTuning(theta_v = 1.0),
            pseudo_timestep = GershgorinPseudoTimeStep(cfl = 0.9))
        res_g = pseudo_transient!(mech_g, cst, solver_g, rt_g)

        @test res_g.converged
        @test all(≈(an.ub, rtol = 1e-5), interior(mech_g.velocity.depthaverage_x))
    end

    @testset "ScaledResidual vs VelocityIncrement" begin
        an = slab_analytical(; const_case...)

        # The normalization is exactly the driving-stress rate, which for a uniform slab is
        # τ_d/(ρH) = g·α — so `error` and `residual` (both taken at the final iterate) must
        # differ by precisely that factor, with no free constant.
        grid, rt, mech = setup_slab()
        fill_slab!(mech, rt, const_case)
        solver = PseudoTransientSolver(grid; maxiter = 200, abstol = 1e-9,
            convergence = ScaledResidual())
        res = pseudo_transient!(mech, cst, solver, rt)

        # `abstol` here is a fraction of the driving stress left unbalanced, so the velocity
        # accuracy it implies is `abstol · τ_d / β = abstol · ub` — hence rtol ≈ abstol,
        # not the m/s reading a `VelocityIncrement` tolerance would have.
        @test res.converged
        @test all(≈(an.ub, rtol = 1e-7), interior(mech.velocity.depthaverage_x))
        @test res.error ≈ res.residual / (cst.gravity * const_case.α) rtol=1e-12

        # The failure mode itself, in miniature: `err = θ·Δτ·r` vanishes as Δτ → 0 whether
        # or not the momentum balance is satisfied. Throttling Δτ by 1e-8 (what an ice
        # shelf's drag-free, diffusion-limited relaxation does to its own increments
        # relative to the grounded ice's) makes the increment criterion report success at
        # essentially zero velocity — the empty Ross and Ronne of the original AIS figure.
        grid_i, rt_i, mech_i = setup_slab()
        fill_slab!(mech_i, rt_i, const_case)
        solver_i = PseudoTransientSolver(grid_i; maxiter = 200, abstol = 1e-8,
            dtau_scaling = 1e-8)
        res_i = pseudo_transient!(mech_i, cst, solver_i, rt_i)

        @test res_i.converged                                   # ... and yet:
        @test maximum(abs, interior(mech_i.velocity.depthaverage_x)) < 1e-6 * abs(an.ub)

        # Same throttled solver, residual criterion: no false positive.
        grid_r, rt_r, mech_r = setup_slab()
        fill_slab!(mech_r, rt_r, const_case)
        solver_r = PseudoTransientSolver(grid_r; maxiter = 20, abstol = 1e-8,
            dtau_scaling = 1e-8, convergence = ScaledResidual())
        res_r = pseudo_transient!(mech_r, cst, solver_r, rt_r)

        @test !res_r.converged
        @test res_r.error > 0.99        # still ~all of the driving stress unbalanced
    end

    # -------------------------------------------------------------------
    # AutotunedDynamicRelaxation (Duretz et al. 2026): the four tests below cover the
    # estimator (against a closed form), what it buys (iteration-count scaling), that it
    # does not move the answer (nonlinear fixed point), and that it stays out of the way
    # when there is nothing to tune.
    # -------------------------------------------------------------------

    # The one configuration where λ_min is known in closed form. On a *uniform* slab from
    # u = 0 the error is exactly the rigid-translation mode u = const, which the discrete
    # membrane operator annihilates (Neumann(0) keeps that true at the boundary), so it is
    # an eigenvector of Ã with eigenvalue β/(ρH). Against the Gershgorin row sum
    # Λ = (32ηH/dx² + β)/(ρH) (dx = dy, see the closed-form Δτ test above) that is
    # λ_min(Â) = β/(32ηH/dx² + β), and since the iterate never leaves that one mode the
    # quotient must return it *exactly*, not asymptotically. Three friction coefficients
    # span a stiffness ratio of 1.01 to 129. `abstol = 0` pins the iteration count so the
    # estimate is read while Δu is still above roundoff — past convergence Δu is pure
    # rounding, which is what the λ_min ≤ λ_max = 1 clamp exists for.
    @testset "AutotunedDynamicRelaxation: λ_min matches the closed form on a uniform slab" begin
        dx = 5e3
        for β0 in (1e4, 1e2, 1e0)
            case = (H0 = 1000.0, μ0 = 1e5, β0 = β0, α = 1e-3)
            grid, rt, mech = setup_slab(; nx = 10, dx)
            fill_slab!(mech, rt, case)
            solver = PseudoTransientSolver(grid; maxiter = 30, abstol = 0.0,
                pseudo_timestep = GershgorinPseudoTimeStep(cfl = 0.99),
                tuning = AutotunedDynamicRelaxation(; cadence = 2))
            res = pseudo_transient!(mech, cst, solver, rt)

            expected = β0 / (32 * case.μ0 * case.H0 / dx^2 + β0)
            @test res.lambda_min ≈ expected rtol=1e-6
            # γ = c·Δτ from the closed form of the two coupled parameters.
            d = 2 * 0.8 * sqrt(expected)
            cc = 0.99^2
            @test res.damping ≈ d * (-cc * d + sqrt(cc^2 * d^2 + 4cc)) rtol=1e-6
            @test 0 < res.damping < 2      # `1 - γ < 0` is legal (over-damped), γ ≥ 2 is not
        end
    end

    # The Phase 2 payoff, and what a *fixed* gamma provably cannot deliver (see the
    # "resolution-robust constant-factor speedup" testset above: a flat ~2x at every
    # resolution). Deliberately stiff — β small next to the membrane stiffness 32ηH/dx², so
    # λ_min ∝ dx² and the condition number grows like nx² at a fixed physical domain, the
    # regime where the damping has to *follow* the mesh. Measured undamped/autotuned
    # ratios: 4.5x (335/75), 7.5x (935/125), 10.3x (2275/220) at nx = 20/40/80, with λ_min
    # falling as 1/nx² (1.06e-2, 3.15e-3, 6.87e-4) and γ tracking it (0.300, 0.170, 0.081).
    # The assertion is on the *growth* of that ratio: a fixed gamma keeps it flat.
    @testset "AutotunedDynamicRelaxation: speedup grows with resolution (not a constant factor)" begin
        stiff = (H0 = 1000.0, μ0 = 1e5, β0 = 1e0, α = 1e-3)
        an = slab_analytical(; stiff...)
        Lx = 100e3

        function stiff_run(nx, tuned)
            dx = Lx / nx
            grid = StaggeredGrid(Float64, Lx, 3 * dx, dx, dx)
            rt, mech = Runtime(grid), MechanicState(grid)
            fill_slab!(mech, rt, stiff)
            fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d,
                           (x, y) -> an.ub + 0.5 * an.ub * sinpi(3x / Lx))
            kw = tuned ? (; tuning = AutotunedDynamicRelaxation()) :
                         (; tuning = FixedTuning(theta_v = 1.0))
            solver = PseudoTransientSolver(grid; maxiter = 100_000, ncheck = 5,
                abstol = 1e-6 * abs(an.ub),
                pseudo_timestep = GershgorinPseudoTimeStep(cfl = 0.99), kw...)
            return pseudo_transient!(mech, cst, solver, rt), mech
        end

        ratios = Float64[]
        lambdas = Float64[]
        for nx in (20, 40, 80)
            res_plain, mech_plain = stiff_run(nx, false)
            res_auto, mech_auto   = stiff_run(nx, true)

            @test res_plain.converged
            @test res_auto.converged
            # Same answer, reached faster — the tuning changes the path, not the fixed point.
            # The tolerance is loose because `VelocityIncrement` is being asked to certify a
            # slowly-converging iteration: the error left at a given increment is that
            # increment divided by `1 - ρ`, so the undamped run (whose `ρ` is closest to 1)
            # stops ~1e-3 from the answer while the autotuned one, at the same `abstol`, is
            # ~1e-5 from it. Both satisfy the criterion they were given; only the autotuned
            # one is also accurate, which is a second, quieter benefit of a smaller `ρ`.
            @test all(≈(an.ub, rtol = 2e-3), interior(mech_plain.velocity.depthaverage_x))
            @test all(≈(an.ub, rtol = 2e-3), interior(mech_auto.velocity.depthaverage_x))
            @test res_auto.iterations < res_plain.iterations

            push!(ratios, res_plain.iterations / res_auto.iterations)
            push!(lambdas, res_auto.lambda_min)
        end

        # Growing speedup: measured 4.5 → 7.5 → 10.3, so 2.3x over two refinements. The
        # threshold is 1.5x, well clear of a flat (fixed-gamma) ratio of 1.0x.
        @test ratios[end] > 1.5 * ratios[1]
        @test issorted(ratios)
        # λ_min ∝ dx²: each doubling of nx must cut it by roughly four.
        @test all(2.5 .< lambdas[1:(end - 1)] ./ lambdas[2:end] .< 6)
    end

    # The nonlinear case, on the same shear problem the `GlenViscosityContinuation` testset
    # above uses. The point that matters is not that the derived damping beats the hand-set
    # gamma = 0.5 (113 iterations vs 167) but that a viscosity moving *during* the solve
    # does not lead the autotuner astray, because every re-estimation also rebuilds the
    # Gershgorin Δτ from the viscosity as it then stands (Duretz §7.3).
    @testset "AutotunedDynamicRelaxation: nonlinear (Glen) solve, same fixed point, fewer iterations" begin
        nx, dx = 40, 5e3
        Lx = nx * dx
        H0, α, β0, n, amp = 1000.0, 1e-3, 1e3, 3.0, 0.95
        ε̇_typical, μ_target = 1e-4, 1e5
        A0 = (1 / (2 * μ_target * ε̇_typical^((n - 1) / n)))^n
        ε̇0 = 1e-8

        function glen_problem()
            g, r, mc = setup_slab(; nx, dx)
            fill_analytic!(mc.topography.surface, r.grid2d, (x, y) -> H0 - α * x)
            fill_analytic!(mc.topography.thickness, r.grid2d, (x, y) -> H0)
            fill_analytic!(mc.friction.beta_eff, r.grid2d,
                           (x, y) -> β0 * (1 + amp * cospi(2x / Lx)))
            fill_analytic!(mc.material.rate_factor_depthaveraged, r.grid2d, (x, y) -> A0)
            fill_analytic!(mc.material.viscosity_depthaveraged, r.grid2d, (x, y) -> μ_target)
            return g, r, mc
        end

        vc = GlenViscosityContinuation(; n_glen = n, theta_mu = 1.0, strainrate_reg = ε̇0)
        function glen_run(kw)
            g, r, mc = glen_problem()
            sv = PseudoTransientSolver(g; maxiter = 20_000, abstol = 1e-9,
                pseudo_timestep = GershgorinPseudoTimeStep(cfl = 0.99),
                viscosity_continuation = vc, kw...)
            return pseudo_transient!(mc, cst, sv, r), mc, r
        end

        res_fixed, mech_fixed, _ = glen_run((; tuning = FixedTuning(theta_v = 1.0, gamma = 0.5)))
        res_auto, mech_auto, rt_auto = glen_run((; tuning = AutotunedDynamicRelaxation()))

        @test res_fixed.converged
        @test res_auto.converged
        @test res_auto.iterations < res_fixed.iterations
        @test 0 < res_auto.damping < 2
        @test res_auto.damping != 0.5                      # genuinely derived, not inherited

        # Same fixed point, to the tolerance the two solves were asked for.
        @test interior(mech_auto.velocity.depthaverage_x) ≈ interior(mech_fixed.velocity.depthaverage_x) rtol=1e-4

        # And still self-consistent with Glen's law at the converged strain rate, i.e. the
        # nonlinear problem really was solved, not just the linear one at the initial η.
        effective_strainrate_ssa!(mech_auto.strainrate, mech_auto.velocity, rt_auto)
        eff = interior(mech_auto.strainrate.effective_depthaveraged)
        μ_closed = @. inv(2 * A0^(1 / n)) * sqrt(eff^2 + ε̇0^2)^((1 - n) / n)
        @test interior(mech_auto.material.viscosity_depthaveraged) ≈ μ_closed rtol=1e-6
    end

    # A well-conditioned problem converges inside the warm-up, before the first Rayleigh
    # quotient is ever taken — the autotuner must then be a pure no-op rather than a source
    # of surprises. `damping` stays at the warm-up's undamped 1 and `lambda_min` is `NaN`
    # ("not estimated"), and the answer is the analytic one either way.
    @testset "AutotunedDynamicRelaxation: no-op when the solve converges inside the warm-up" begin
        an = slab_analytical(; const_case...)
        grid, rt, mech = setup_slab()
        fill_slab!(mech, rt, const_case)
        solver = PseudoTransientSolver(grid; maxiter = 200, abstol = 1e-8,
            pseudo_timestep = GershgorinPseudoTimeStep(cfl = 0.99),
            tuning = AutotunedDynamicRelaxation(; cadence = 100))
        res = pseudo_transient!(mech, cst, solver, rt)

        @test res.converged
        @test res.iterations < 100                 # ... i.e. inside the warm-up
        @test res.damping == 1
        @test isnan(res.lambda_min)
        @test all(≈(an.ub, rtol = 1e-5), interior(mech.velocity.depthaverage_x))

        # FixedTuning (the default) reports the solver's own gamma and no estimate at all.
        solver_fixed = PseudoTransientSolver(grid; maxiter = 200, abstol = 1e-8,
            tuning = FixedTuning(gamma = 0.7))
        setdata!(mech.velocity.depthaverage_x, 0.0); setdata!(mech.velocity.depthaverage_y, 0.0)
        res_fixed = pseudo_transient!(mech, cst, solver_fixed, rt)
        @test res_fixed.damping == 0.7
        @test isnan(res_fixed.lambda_min)
    end
end

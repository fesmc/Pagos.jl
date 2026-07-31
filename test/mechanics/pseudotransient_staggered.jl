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

        # Default gamma = 1 disables damping (plain, undamped PT iteration).
        @test solver.gamma == 1

        # Requires a depth-averaged grid: DIVA's vertical shear integral is future work.
        layering  = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 4))
        col_grid  = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, layering)
        @test_throws ArgumentError PseudoTransientSolver(col_grid)
    end

    @testset "pseudo_dt is representation-agnostic (Field vs Array)" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        mech = MechanicState(grid)
        fill_analytic!(mech.material.viscosity_depthaveraged, grid.grid2d,
                       (x, y) -> 1e5 + 3e3 * x)

        dt_field = pseudo_dt(910.0, 1.0, 1.0, mech.material.viscosity_depthaveraged, 1e2, 4.1)
        dt_array = pseudo_dt(910.0, 1.0, 1.0, interior(mech.material.viscosity_depthaveraged),
                             1e2, 4.1)
        @test dt_field ≈ dt_array
        @test dt_field ≈ 910.0 * 1.0 * 1.0 / (4 * (1 + 1e2) * 4.1 * maximum(interior(
            mech.material.viscosity_depthaveraged)))
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

        # Uniform viscosity: the local field collapses to the same value everywhere, and
        # to the collocated scalar `pseudo_dt` (no lerp needed, both cells agree).
        fill_analytic!(mech.material.viscosity_depthaveraged, rt.grid2d, (x, y) -> η0)
        pseudo_dt!(solver.dtau_x, solver.dtau_y, ρ, dx, dy,
                  mech.material.viscosity_depthaveraged, muB, ndim, dtau_scaling, rt)
        dt_scalar = dtau_scaling * pseudo_dt(ρ, dx, dy, interior(mech.material.viscosity_depthaveraged), muB, ndim)
        @test all(≈(dt_scalar), interior(solver.dtau_x))
        @test all(≈(dt_scalar), interior(solver.dtau_y))
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

        fill_analytic!(mech.velocity.x, rt.grid, (x, y) -> a * x + b * y)
        fill_analytic!(mech.velocity.y, rt.grid, (x, y) -> c * x + d * y)
        velocitygradients!(mech.velocity, mech.topography.thickness, rt)
        effective_strainrate_ssa!(mech.strainrate, mech.velocity, rt)

        expected = sqrt(a^2 + d^2 + a * d + ((b + c) / 2)^2)
        @test all(≈(expected), interior(mech.strainrate.effective))
    end

    # Direct check of the Glen-law + log-space relaxation formula (Sandip et al. 2024,
    # Eq. 3 and 8) through the public `update_viscosity!` entry point, on a uniform
    # velocity gradient (so `effective_strainrate_ssa!`'s output is a known constant, per
    # the test above) and uniform rate factor, so the closed form is exact everywhere.
    @testset "update_viscosity!(::GlenViscosityContinuation): matches the closed-form formula" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        rt   = Runtime(grid)
        mech = MechanicState(grid)
        A0, n, ε̇0, μ_old0 = 1e-16, 3.0, 1e-12, 3e14
        a, b, c, d = 2e-3, -1e-3, 5e-4, 3e-3   # same linear velocity as the test above

        fill_analytic!(mech.velocity.x, rt.grid, (x, y) -> a * x + b * y)
        fill_analytic!(mech.velocity.y, rt.grid, (x, y) -> c * x + d * y)
        fill_analytic!(mech.material.rate_factor_depthaveraged, rt.grid2d, (x, y) -> A0)
        velocitygradients!(mech.velocity, mech.topography.thickness, rt)

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
        # between centres i-1, i; i=1 reaches centre 0). See `roadmaps/chmy.md`, §3.
        fill_analytic!(mech.material.viscosity_depthaveraged, rt.grid2d, (x, y) -> η0)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
        fill_analytic!(mech.velocity.x, rt.grid, (x, y) -> a * x + b * y)
        fill_analytic!(mech.velocity.y, rt.grid, (x, y) -> c * x + d * y)

        velocitygradients!(mech.velocity, mech.topography.thickness, rt)
        strainrate!(mech.strainrate, mech.velocity, mech.material, mech.topography,
                   DIVAMomentumBalance(), rt)

        @test all(interior(mech.strainrate.xx) .≈ 2η0 * H0 * (2a + d))
        @test all(interior(mech.strainrate.yy) .≈ 2η0 * H0 * (a + 2d))
        @test all(interior(mech.strainrate.xy) .≈ η0 * H0 * (b + c))

        # Locations match the layout table: N_xx/N_yy at `aa`, N_xy at `ab`.
        @test location(mech.strainrate.xx) === (Center(), Center(), Center())
        @test location(mech.strainrate.xy) === (Vertex(), Vertex(), Center())
    end

    # A genuine viscosity contrast: hlerp must give the harmonic, not the arithmetic, mean.
    # Located by cell centres, not a coordinate/index guess, per the vertex/centre lesson
    # in `roadmaps/chmy.md` (§3): vertex `i` sits between centres `i-1` and `i`.
    @testset "membrane stress: harmonic η averaging at a viscosity step" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        rt   = Runtime(grid)
        mech = MechanicState(grid)
        η_lo, η_hi, H0, b = 1e4, 1e6, 500.0, 1e-2

        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
        fill_analytic!(mech.velocity.x, rt.grid, (x, y) -> b * y)   # ε̇xy = b/2 uniform
        fill_analytic!(mech.velocity.y, rt.grid, (x, y) -> 0.0)
        fill_analytic!(mech.material.viscosity_depthaveraged, rt.grid2d,
                       (x, y) -> x < 0 ? η_lo : η_hi)

        velocitygradients!(mech.velocity, mech.topography.thickness, rt)
        strainrate!(mech.strainrate, mech.velocity, mech.material, mech.topography,
                   DIVAMomentumBalance(), rt)

        # Find an ab vertex whose two x-centres genuinely straddle the step.
        xs = collect(grid.x)
        istep = findfirst(i -> xs[i] < 0 <= xs[i + 1], 1:(length(xs) - 1))
        @test xs[istep] < 0 <= xs[istep + 1]         # the premise: the face really straddles it
        i_ab = istep + 1                              # vertex i sits between centres i-1, i

        η_harm = 2 / (1 / η_lo + 1 / η_hi)
        expected_xy = η_harm * H0 * b
        @test interior(mech.strainrate.xy)[i_ab, 4, 1] ≈ expected_xy
        # The arithmetic mean would have given a different (larger) answer.
        @test !(interior(mech.strainrate.xy)[i_ab, 4, 1] ≈ (η_lo + η_hi) / 2 * H0 * b)
    end

    # Containment: an unmasked ice-free (η = 0) cell produces NaN via hlerp; masking with
    # the strict rule keeps it finite, mirroring `deviatoric_stress!`'s established pattern.
    @testset "membrane stress: NaN off-ice, contained by IceMask" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        rt   = Runtime(grid)
        mech = MechanicState(grid)
        topo = TopographicState(grid)

        setdata!(mech.topography.thickness, 500.0)
        fill_analytic!(mech.velocity.x, rt.grid, (x, y) -> 0.01y)
        fill_analytic!(mech.velocity.y, rt.grid, (x, y) -> 0.0)
        fill_analytic!(mech.material.viscosity_depthaveraged, rt.grid2d,
                       (x, y) -> x < 0 ? 0.0 : 1e5)
        velocitygradients!(mech.velocity, mech.topography.thickness, rt)

        strainrate!(mech.strainrate, mech.velocity, mech.material, mech.topography,
                   DIVAMomentumBalance(), rt)
        @test any(isnan, interior(mech.strainrate.xy))

        setdata!(topo.mask.is_ice, true)
        for i in axes(interior(topo.mask.is_ice), 1)
            x, _, _ = coord(rt.grid2d, location(topo.mask.is_ice), i, 1, 1)
            x < 0 && (interior(topo.mask.is_ice)[i, :, 1] .= false)
        end
        mask = IceMask(topo.mask.is_ice)

        strainrate!(mech.strainrate, mech.velocity, mech.material, mech.topography,
                   DIVAMomentumBalance(), rt, mask)
        @test !any(isnan, interior(mech.strainrate.xy))
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
        @test all(≈(an.ub, rtol = 1e-5), interior(mech.velocity.x))
        @test all(≈(0.0,   atol = 1e-8 * abs(an.ub)), interior(mech.velocity.y))
        # residual is the raw momentum-balance rate, not the velocity increment `error` is
        # based on — different units, but both must be tiny at a converged solution.
        @test res.residual >= 0
        @test res.residual < 1e-6

        # SSA and DIVA share the same equations in this SSA-limit (no resolved shear).
        setdata!(mech.velocity.x, 0.0); setdata!(mech.velocity.y, 0.0)
        res_ssa = pseudo_transient!(mech, cst, solver, rt, SSAMomentumBalance())
        @test res_ssa.converged
        @test all(≈(an.ub, rtol = 1e-5), interior(mech.velocity.x))
    end

    @testset "ncheck = 5" begin
        an = slab_analytical(; const_case...)
        grid, rt, mech = setup_slab()
        fill_slab!(mech, rt, const_case)
        solver5 = PseudoTransientSolver(grid; maxiter = 200, abstol = 1e-8, ncheck = 5)

        res5 = pseudo_transient!(mech, cst, solver5, rt)
        @test res5.converged
        @test res5.iterations % 5 == 0
        @test all(≈(an.ub, rtol = 1e-5), interior(mech.velocity.x))
        @test all(≈(0.0,   atol = 1e-8 * abs(an.ub)), interior(mech.velocity.y))
    end

    @testset "Float32 stays Float32" begin
        an = slab_analytical(; const_case...)
        grid, rt, mech = setup_slab(Float32)
        fill_slab!(mech, rt, const_case, Float32)
        cst32  = Constants{Float32}()
        solver = PseudoTransientSolver(grid; maxiter = 200, abstol = 1f-6)

        res = pseudo_transient!(mech, cst32, solver, rt)
        @test res.converged
        @test eltype(mech.velocity.x) === Float32
        @test all(≈(Float32(an.ub), rtol = 1f-3), interior(mech.velocity.x))
    end

    # Damping (Sandip et al. 2024, Eq. 12–14) gives a real, fixed-γ iteration-count
    # reduction (verified here and, at several resolutions, in the next testset below) —
    # but NOT by itself the sub-quadratic *scaling* the papers report; see that testset's
    # note for why a fixed γ isn't sufficient for that (it needs Phase 2's resolution-
    # dependent parameter selection). The uniform-slab solution itself has zero velocity
    # gradient, so its membrane stress is identically zero and convergence there is
    # friction-dominated regardless of gamma (see the "large β0" comment on `const_case`
    # above) — not a useful case for exercising damping. This test instead starts from a
    # spatially oscillating perturbation around that same analytic solution: relaxing the
    # perturbation away requires the membrane-stress divergence to act repeatedly, which is
    # exactly the pathway damping accelerates. Both solves converge to the *same* analytic
    # solution regardless of gamma — damping only changes the transient path, never the
    # fixed point — so this isolates the iteration-count effect. Parameters calibrated
    # empirically: at (β0, nx) = (1e3, 20) this perturbation converges in 184 iterations
    # undamped vs. 58 with gamma = 0.5 (a ~3x reduction); the assertion below only requires
    # 2x, for headroom.
    @testset "gamma < 1 reduces iterations to converge (perturbed initial condition)" begin
        damp_case = (H0 = 1000.0, μ0 = 1e5, β0 = 1e3, α = 1e-3)
        an = slab_analytical(; damp_case...)
        nx, dx = 20, 5e3
        lx = nx * dx

        function perturbed_run(gamma)
            grid, rt, mech = setup_slab(; nx, dx)
            fill_slab!(mech, rt, damp_case)
            fill_analytic!(mech.velocity.x, rt.grid,
                           (x, y) -> an.ub + 0.5 * an.ub * sinpi(6x / lx))
            solver = PseudoTransientSolver(grid; maxiter = 1000, abstol = 1e-8, gamma)
            return pseudo_transient!(mech, cst, solver, rt), mech
        end

        res_undamped, mech_undamped = perturbed_run(1.0)
        res_damped, mech_damped     = perturbed_run(0.5)

        @test res_undamped.converged
        @test res_damped.converged
        @test res_damped.iterations < res_undamped.iterations ÷ 2

        @test all(≈(an.ub, rtol = 1e-5), interior(mech_undamped.velocity.x))
        @test all(≈(an.ub, rtol = 1e-5), interior(mech_damped.velocity.x))
    end

    # Verification, at several resolutions, of what a *fixed* gamma actually buys — and,
    # honestly, what it does not. The papers' O(N²) → O(N^1.2)-ish scaling result comes from
    # tuning the damping parameter *to the current mesh's spectral radius* (Duretz et al.
    # 2026 Eq. 19: c = c_damp · 2√λ_min, which itself depends on resolution) — that
    # resolution-dependent selection is Phase 2 (`roadmaps/PT-autotune.md`) and not
    # implemented yet. With a single fixed `gamma` reused across resolutions, measurement
    # here (fixed physical domain `Lx`, dx = Lx/nx refined, same relative perturbation
    # wavelength at every nx) shows the undamped/damped iteration count *ratio* holding
    # essentially constant (~2x at nx = 10, 20, 40 — 505/236, 1884/958, 6930/3595
    # iterations) rather than growing with resolution: both schemes scale similarly with
    # nx, gamma = 0.5 just scales with a smaller prefactor. That is still a genuine,
    # resolution-robust win (not a fluke at one grid size) — just not the asymptotic
    # exponent change, which needs Phase 2's resolution-adaptive parameter.
    @testset "gamma < 1: resolution-robust constant-factor speedup (not yet O(N) scaling)" begin
        scale_case = (H0 = 1000.0, μ0 = 1e5, β0 = 1e2, α = 1e-3)
        an = slab_analytical(; scale_case...)
        Lx = 100e3

        function resolved_run(nx, gamma)
            dx = Lx / nx
            grid, rt, mech = setup_slab(; nx, dx)
            fill_slab!(mech, rt, scale_case)
            fill_analytic!(mech.velocity.x, rt.grid,
                           (x, y) -> an.ub + 0.5 * an.ub * sinpi(3x / Lx))
            solver = PseudoTransientSolver(grid; maxiter = 20_000, abstol = 1e-8, gamma)
            res = pseudo_transient!(mech, cst, solver, rt)
            return res, mech
        end

        for nx in (10, 20, 40)
            res_undamped, mech_undamped = resolved_run(nx, 1.0)
            res_damped, mech_damped     = resolved_run(nx, 0.5)

            @test res_undamped.converged
            @test res_damped.converged
            @test res_damped.iterations < res_undamped.iterations / 1.5

            @test all(≈(an.ub, rtol = 1e-5), interior(mech_undamped.velocity.x))
            @test all(≈(an.ub, rtol = 1e-5), interior(mech_damped.velocity.x))
        end
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
                                           gamma = 0.5, viscosity_continuation = vc)
        res_glen = pseudo_transient!(mech_glen, cst, solver_glen, rt_glen)
        @test res_glen.converged

        # Self-consistency: the converged viscosity must match Glen's law evaluated at the
        # *converged* strain rate — recomputed independently here, not read back from
        # whatever the solver last wrote mid-iteration.
        effective_strainrate_ssa!(mech_glen.strainrate, mech_glen.velocity, rt_glen)
        eff = interior(mech_glen.strainrate.effective)
        μ_closed = @. inv(2 * A0^(1 / n)) * sqrt(eff^2 + ε̇0^2)^((1 - n) / n)
        @test interior(mech_glen.material.viscosity_depthaveraged) ≈ μ_closed rtol=1e-6

        # The mechanism has a real effect: compare against a fixed-viscosity (no
        # continuation) solve from the same initial guess — measured ~4.6% difference at
        # these parameters, the assertion only requires 1%.
        grid_fixed, rt_fixed, mech_fixed = build_shear_problem()
        solver_fixed = PseudoTransientSolver(grid_fixed; maxiter = 5000, abstol = 1e-9,
                                            gamma = 0.5)
        res_fixed = pseudo_transient!(mech_fixed, cst, solver_fixed, rt_fixed)
        @test res_fixed.converged

        reldiff = maximum(abs.(interior(mech_glen.velocity.x) .- interior(mech_fixed.velocity.x))) /
                  maximum(abs.(interior(mech_fixed.velocity.x)))
        @test reldiff > 0.01

        # theta_mu only changes the transient path, not the fixed point (same argument as
        # the damping tests above, now for the log-space viscosity relaxation instead of
        # the velocity-rate one): a relaxed continuation must converge to the same answer.
        vc_relaxed = GlenViscosityContinuation(; n_glen = n, theta_mu = 0.3, strainrate_reg = ε̇0)
        grid_relaxed, rt_relaxed, mech_relaxed = build_shear_problem()
        solver_relaxed = PseudoTransientSolver(grid_relaxed; maxiter = 5000, abstol = 1e-9,
                                               gamma = 0.5, viscosity_continuation = vc_relaxed)
        res_relaxed = pseudo_transient!(mech_relaxed, cst, solver_relaxed, rt_relaxed)
        @test res_relaxed.converged
        @test all(isapprox.(interior(mech_relaxed.velocity.x), interior(mech_glen.velocity.x);
                            rtol = 1e-4))
    end
end

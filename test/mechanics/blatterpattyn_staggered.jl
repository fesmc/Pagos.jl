using Pagos
using Test

include("../test_helpers/chmy.jl")

# Phase 1 of `roadmaps/blatter-pattyn.md`: the explicit Blatter-Pattyn pseudo-transient
# solver, reusing the SSA/DIVA iteration skeleton on the column grid. Verified the same way
# `pseudotransient_staggered.jl` verifies SSA/DIVA: closed-form pieces first, then the full
# `pseudo_transient!` solve against the one case with an analytic solution (§Phase 1 of the
# roadmap) — never against the 2D solve, which discretizes different points entirely.

# -----------------------------------------------------------------------
# Analytical uniform-slab solution, constant (non-Glen) viscosity.
#
# A slab uniform in x/y collapses the membrane terms (∂x, ∂y of a spatially uniform field
# vanish), leaving the pure vertical-shear ODE `∂z(µ ∂u/∂z) = ρg ∂s/∂x` with `µ` constant,
# bed flux `µ ∂u/∂z|_b = β u_b` and stress-free surface `µ ∂u/∂z|_s = 0`. Integrating twice
# (surface BC fixes the integration constant, bed BC fixes u_b) gives, with `ζ = z/H`:
#
#   u(ζ) = u_b + (ρ g α H²/µ) · (ζ − ζ²/2),     u_b = ρ g H α / β
#
# `u_b` and the depth average `ū = u_b + ρgH²α/(3µ)` are exactly `slab_analytical`'s SSA
# closed form (`test/mechanics/pseudotransient_staggered.jl`) — expected, since Robinson et
# al. derive SSA *from* BP by depth-integrating exactly this balance, so BP's z-integral
# must reproduce it on a geometry where SSA itself is exact. This is the one case with a
# closed form: it exercises the vertical operator and both boundary conditions, and nothing
# else (`roadmaps/blatter-pattyn.md`, Phase 1 verification).
# -----------------------------------------------------------------------

function bp_slab_analytical(; H0, μ0, β0, α, ρ = 910.0, g = 9.81)
    ub = ρ * g * H0 * α / β0
    ubar = ub + ρ * g * H0^2 * α / (3μ0)
    profile(ζ) = ub + (ρ * g * α * H0^2 / μ0) * (ζ - ζ^2 / 2)
    return (; ub, ubar, profile)
end

# How wide the slab domain has to be for the closed form to be the answer in the middle.
#
# The membrane terms vanish in the *continuous* slab problem, but not in the discrete one:
# the domain edges carry the `Neumann(0)` halo placeholder (`roadmaps/chmy.md` Phase 4)
# rather than a real per-equation boundary condition, so a small membrane residual is seeded
# at x = 0 and x = L. This used to be invisible because `σxx` was identically zero on a
# slab; with the terrain-following metric correction in `velocitygradients!`
# (`terrain_metric_correction!`, `roadmaps/blatter-pattyn-equations.md` A1) the surface slope
# makes `σxx = 4µ α ∂u/∂z ≠ 0`, and the edge artifact becomes visible.
#
# It decays away from the edges on the membrane (SSA coupling) length scale
#
#   L_m = sqrt(4 µ H / β) = sqrt(4 · 1e8 · 1e3 / 1e3) = 20 km = 20 cells at Δx = 1 km,
#
# so the original 8-cell domain was *entirely* boundary layer and the closed form was never
# recovered anywhere in it. At 128 cells the midpoint sits 3.2 L_m from either edge and the
# midcolumn error drops to 9.7e-5 — an order of magnitude under the 1e-3 tolerance, and at
# the ~8e-5 floor set by the vertical discretization itself (measured: 1.1e-3 at nx = 8,
# 1.7e-4 at 64, 9.7e-5 at 128, 7.9e-5 at 256). Do not shrink this back without also giving
# the momentum solver real domain-edge boundary conditions.
const SLAB_NX = 128

@testset "Blatter-Pattyn pseudo-transient momentum solver (C-grid staggered)" begin
    cst = Constants{Float64}()
    layering(T = Float64; nz = 4) = CorrectedVerticalLayering(T, QuadraticSigmaTransform(T, nz))

    @testset "PseudoTransientSolver(::StaggeredGrid, ::BlatterPattynMomentumBalance)" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, layering())
        momentum = BlatterPattynMomentumBalance()
        solver = PseudoTransientSolver(grid, momentum)

        # Work arrays are ACX3/ACY3 on `grid.grid`, matching `mech.velocity.x`/`y` — not
        # `grid.grid2d`, unlike the depth-averaged constructor.
        @test location(solver.velocity_x_old) === (Vertex(), Center(), Center())  # acx
        @test location(solver.velocity_y_old) === (Center(), Vertex(), Center())  # acy
        @test size(interior(solver.velocity_x_old)) == (grid.nx + 1, grid.ny, grid.nz)
        @test size(interior(solver.velocity_y_old)) == (grid.nx, grid.ny + 1, grid.nz)
        @test size(interior(solver.dtau_x)) == size(interior(solver.velocity_x_old))
        @test size(interior(solver.dtau_y)) == size(interior(solver.velocity_y_old))
        @test eltype(solver.velocity_x_old) === Float64

        # The depth-averaged constructor is unaffected by this method's existence — same
        # call, same shapes as it always had.
        solver2d = PseudoTransientSolver(grid)
        @test size(interior(solver2d.velocity_x_old)) == (grid.nx + 1, grid.ny, 1)
    end

    @testset "momentum balance grid guards" begin
        flat = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        col  = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, layering())
        momentum = BlatterPattynMomentumBalance()

        mech_flat, rt_flat = MechanicState(flat), Runtime(flat)
        solver = PseudoTransientSolver(flat, momentum; maxiter = 1)
        @test_throws ArgumentError pseudo_transient!(mech_flat, cst, solver, rt_flat,
                                                      momentum)

        mech_col, rt_col = MechanicState(col), Runtime(col)
        solver_col = PseudoTransientSolver(col, momentum; maxiter = 1)
        # Doesn't throw — reaches the (short) iteration loop instead.
        res = pseudo_transient!(mech_col, cst, solver_col, rt_col, momentum)
        @test res.iterations == 1
    end

    @testset "effective_strainrate_bp!: exact under a linear velocity+shear field" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, layering())
        rt = Runtime(grid)
        mech = MechanicState(grid)

        a, b, c, d = 0.1, 0.2, -0.05, 0.3
        e, f = 0.02, -0.01
        H0 = 100.0
        fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> a * x + b * y + e * ζ)
        fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> c * x + d * y + f * ζ)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)

        velocitygradients!(mech.velocity, mech.topography.thickness, rt)
        effective_strainrate_bp!(mech, rt)

        uz, vz = e / H0, f / H0
        expected = sqrt(a^2 + d^2 + a * d + ((b + c) / 2)^2 + uz^2 / 4 + vz^2 / 4)
        @test all(≈(expected, rtol = 1e-10), interior(mech.strainrate.effective))
    end

    @testset "membranestress!(::MomentumBalance3D): exact under uniform µ" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, layering())
        rt = Runtime(grid)
        mech = MechanicState(grid)

        a, b, c, d = 0.1, 0.2, -0.05, 0.3
        e = 0.02
        H0, μ0 = 100.0, 1e7
        fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> a * x + b * y + e * ζ)
        fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> c * x + d * y)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0)

        velocitygradients!(mech.velocity, mech.topography.thickness, rt)
        membranestress!(mech, BlatterPattynMomentumBalance(), rt)

        @test all(≈(2μ0 * (2a + d)), interior(mech.stress.xx))
        @test all(≈(2μ0 * (a + 2d)), interior(mech.stress.yy))
        @test all(≈(μ0 * (b + c)), interior(mech.stress.xy))
        @test all(≈(μ0 * (e / H0)), interior(mech.stress.xz))
        @test all(≈(0.0, atol = 1e-8), interior(mech.stress.yz))  # v has no z-dependence
    end

    @testset "margin faces keep their vertical operator" begin
        # The regression the AIS 8 km run exposed (`roadmaps/blatter-pattyn.md`, Phase 1
        # "margin σxz"). `NODE_ACX_AC` and `NODE_ACX` resolve to the *same* cell pair, so
        # gating σxz on `node_fully_active` while the unknown is created under `node_active`
        # zeroed the whole vertical operator along the margin — and BP's only tie between a
        # column and the bed is the k = 1 flux, so layers above drift without a restoring
        # force. Half the domain is iced, giving one clean margin at the `acx` face i = m + 1.
        # A realistic aspect ratio (Δx = 1 km against H = 100 m) on purpose: it is what makes
        # the vertical term dominate the Gershgorin row sum, and hence what lets the Δτ
        # assertion below discriminate. At Δx ~ H the horizontal term dominates and a missing
        # vertical operator barely moves Δτ at all.
        nx, ny, nz = 8, 4, 6
        H0, μ0, e, dx = 100.0, 1e7, 0.02, 1e3
        m = nx ÷ 2
        grid = StaggeredGrid(Float64, nx * dx, ny * dx, dx, dx, layering(; nz))
        rt = Runtime(grid)
        mech = MechanicState(grid)
        topo = TopographicState(grid)

        setdata!(topo.mask.is_ice, false)
        interior(topo.mask.is_ice)[1:m, :, 1] .= true
        mask = IceMask(topo.mask.is_ice)

        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
        fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> e * ζ)
        # µ = 0 off ice: the case the strict rule was reached for. A one-sided harmonic mean
        # never inverts it, so σxz must come out finite *and* non-zero on the margin face.
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> 0.0)
        for k in 1:nz, j in -1:(ny + 2), i in 1:m
            mech.material.viscosity[i, j, k] = μ0
        end

        velocitygradients!(mech.velocity, mech.topography.thickness, rt, mask)
        membranestress!(mech, BlatterPattynMomentumBalance(), rt, mask)

        # Interior face (both cells iced): unchanged, the two-way `hlerp`.
        @test mech.stress.xz[m, 2, 3] ≈ μ0 * (e / H0)
        # Margin face i = m + 1 (cell m iced, cell m+1 not): the unknown lives here
        # (`node_active(NODE_ACX)`), so σxz must too — one-sided onto column m.
        @test node_active(mask, Pagos.NODE_ACX, m + 1, 2)
        @test isfinite(mech.stress.xz[m + 1, 2, 3])
        @test mech.stress.xz[m + 1, 2, 3] ≈ μ0 * (e / H0)
        # Fully off-ice face: still zero, and still finite.
        @test mech.stress.xz[m + 3, 2, 3] == 0.0

        # And the Gershgorin bound must see the same µ, or it stops describing the operator
        # it bounds. Both faces are then vertical-stiffness-limited by the *same* µ, so their
        # Δτ agree to within the horizontal row sum — which does legitimately halve at the
        # margin (the ice-free `aa` cell contributes µ = 0), but is ~4 orders of magnitude
        # below the vertical term at this aspect ratio. Drop the vertical term at the margin
        # only, and this ratio jumps from ~1 to ~10⁴.
        solver = PseudoTransientSolver(grid, BlatterPattynMomentumBalance();
                                       tuning = AutotunedDynamicRelaxation())
        fill_analytic!(mech.friction.beta_eff, rt.grid2d, (x, y) -> 1e3)
        pseudo_dt!(solver, mech, cst, rt, BlatterPattynMomentumBalance(), mask)
        @test solver.dtau_x[m + 1, 2, 3] ≈ solver.dtau_x[m, 2, 3] rtol = 1e-2
    end

    @testset "uniform slab, constant viscosity — pseudo_transient!" begin
        const_case = (H0 = 1000.0, μ0 = 1e8, β0 = 1e3, α = 1e-2)
        an = bp_slab_analytical(; const_case...)

        dx = 1e3
        nx = SLAB_NX
        nz = 8
        grid = StaggeredGrid(Float64, nx * dx, 3 * dx, dx, dx, layering(; nz))
        rt = Runtime(grid)
        mech = MechanicState(grid)

        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> const_case.H0)
        fill_analytic!(mech.topography.surface, rt.grid2d, (x, y) -> -const_case.α * x)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> const_case.μ0)
        fill_analytic!(mech.friction.beta_eff, rt.grid2d, (x, y) -> const_case.β0)

        momentum = BlatterPattynMomentumBalance()
        solver = PseudoTransientSolver(grid, momentum;
            maxiter = 20_000, abstol = 1e-11, ncheck = 20,
            tuning = AutotunedDynamicRelaxation())

        res = pseudo_transient!(mech, cst, solver, rt, momentum)
        @test res.converged
        @test res.residual < 1e-6

        # Pointwise column profile against the closed form, at the column furthest from both
        # domain edges — see `SLAB_NX` for why that distance matters.
        i, j = nx ÷ 2, 2
        for k in 1:nz
            ζ = zcenter(rt.grid, k)
            @test mech.velocity.x[i, j, k] ≈ an.profile(ζ) rtol=1e-3
        end
        @test all(≈(0.0, atol = 1e-6 * an.ub), interior(mech.velocity.y))

        # Bed value: the boundary condition folded into `dotvel!`'s vertical-flux assembly.
        @test mech.velocity.x[i, j, 1] ≈ an.profile(zcenter(rt.grid, 1)) rtol=1e-3
    end

    @testset "ScaledResidual convergence: same fixed point as VelocityIncrement" begin
        const_case = (H0 = 1000.0, μ0 = 1e8, β0 = 1e3, α = 1e-2)
        an = bp_slab_analytical(; const_case...)
        dx = 1e3
        grid = StaggeredGrid(Float64, SLAB_NX * dx, 3dx, dx, dx, layering(; nz = 8))
        rt = Runtime(grid)
        mech = MechanicState(grid)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> const_case.H0)
        fill_analytic!(mech.topography.surface, rt.grid2d, (x, y) -> -const_case.α * x)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> const_case.μ0)
        fill_analytic!(mech.friction.beta_eff, rt.grid2d, (x, y) -> const_case.β0)

        momentum = BlatterPattynMomentumBalance()
        solver = PseudoTransientSolver(grid, momentum;
            maxiter = 20_000, abstol = 1e-8, ncheck = 20,
            tuning = AutotunedDynamicRelaxation(), convergence = ScaledResidual())

        res = pseudo_transient!(mech, cst, solver, rt, momentum)
        @test res.converged
        i, j = SLAB_NX ÷ 2, 2
        @test mech.velocity.x[i, j, 1] ≈ an.profile(zcenter(rt.grid, 1)) rtol=1e-3
    end

    @testset "Gershgorin bound reduces to the 2D bound at nz == 1 magnitude order" begin
        # `_check_momentum_grid` rejects `nz == 1` for BP outright, so the reduction the
        # roadmap's §2.1 sanity check calls for cannot be exercised end-to-end through
        # `pseudo_transient!`; this checks the row-sum kernel directly instead, at the
        # smallest legal `nz = 2`, against the closed-form Sandip Eq. 7 scale it must be of
        # the same order as for a well-resolved column (`Δz ≪ Δx`).
        dx = 1e4  # deliberately coarse Δx relative to Δz so the vertical term is not
                  # starved by an unrelated horizontal bound (see the Phase 1 session note
                  # on aspect ratio in `roadmaps/blatter-pattyn.md` §2).
        H0, μ0, β0 = 1000.0, 1e8, 1e3
        grid = StaggeredGrid(Float64, 8dx, 8dx, dx, dx, layering(; nz = 2))
        rt = Runtime(grid)
        mech = MechanicState(grid)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0)
        fill_analytic!(mech.friction.beta_eff, rt.grid2d, (x, y) -> β0)

        momentum = BlatterPattynMomentumBalance()
        solver = PseudoTransientSolver(grid, momentum)
        pseudo_dt!(solver, mech, cst, rt, momentum)

        @test all(>(0), interior(solver.dtau_x))
        @test all(isfinite, interior(solver.dtau_x))
        # Same order of magnitude as the 2D bound at matching µ/H/β/dx (Sandip's slab check
        # in `pseudotransient_staggered.jl`): neither vanishing nor exploding.
        dtau_2d_scale = cst.density_ice * dx^2 / (16μ0)
        @test 1e-4 * dtau_2d_scale < minimum(interior(solver.dtau_x)) < 1e2 * dtau_2d_scale
    end

    @testset "Float32 stays Float32" begin
        const_case = (H0 = 1000.0f0, μ0 = 1f8, β0 = 1f3, α = 1f-2)
        dx = 1f3
        grid = StaggeredGrid(Float32, 8dx, 3dx, dx, dx, layering(Float32))
        rt = Runtime(grid)
        mech = MechanicState(grid)
        cst32 = Constants{Float32}()

        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> const_case.H0)
        fill_analytic!(mech.topography.surface, rt.grid2d, (x, y) -> -const_case.α * x)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> const_case.μ0)
        fill_analytic!(mech.friction.beta_eff, rt.grid2d, (x, y) -> const_case.β0)

        momentum = BlatterPattynMomentumBalance()
        solver = PseudoTransientSolver(grid, momentum;
            maxiter = 20_000, abstol = 1f-5, ncheck = 20,
            tuning = AutotunedDynamicRelaxation())
        res = pseudo_transient!(mech, cst32, solver, rt, momentum)

        @test eltype(mech.velocity.x) === Float32
        @test eltype(solver.dtau_x) === Float32
        @test res.converged
    end
end

using Pagos
using Test

include("../test_helpers/chmy.jl")

# Stage 2 of the DIVA/SSA split (`roadmaps/chmy.md`, Phase 3, decisions 6-12): the pieces
# that give DIVA a genuine vertical-shear correction instead of the SSA-limit-in-disguise
# behaviour Stage 0/1 left it in. Validated against Robinson et al. (2022)'s own equations
# and, for the full solve, their linearised-slab reference state (§3.1) — decision 16, no
# external data.

const _DIVA_PAPER = "Robinson et al. (2022), The Cryosphere 16, 689-709"

@testset "DIVA depth-integrated-viscosity chain (C-grid staggered)" begin

    # ε̇_e² = ūx² + v̄y² + ūxv̄y + ¼(ūy+v̄x)² + ¼(u_z² + v_z²), Eq. 13, with
    # u_z = τ_b,x(1-ζ)/µ(z), Eq. 21 with H cancelled (see the source note on
    # `effective_strainrate_diva!`). Values chosen so the shear term dominates the
    # (small) horizontal one, making the per-layer ζ-dependence visible rather than lost
    # in floating-point noise — a uniform-everything case could pass by accident.
    @testset "effective_strainrate_diva!: matches Eq. 13 pointwise, per layer" begin
        lay  = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 6))
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, lay)
        rt   = Runtime(grid)
        mech = MechanicState(grid)

        a, b, c, d = 0.03, 0.017, -0.011, 0.023
        τbx0, τby0, μ0 = 1e4, -0.6e4, 1e5

        fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d, (x, y) -> a * x + b * y)
        fill_analytic!(mech.velocity.depthaverage_y, rt.grid2d, (x, y) -> c * x + d * y)
        fill_analytic!(mech.stress.base_x, rt.grid2d, (x, y) -> τbx0)
        fill_analytic!(mech.stress.base_y, rt.grid2d, (x, y) -> τby0)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0)
        depthaverage_velocitygradients!(mech.velocity, rt)

        effective_strainrate_diva!(mech, rt)

        horizontal = a^2 + d^2 + a * d + ((b + c) / 2)^2
        for k in 1:grid.nz
            ζ = zcenter(rt.grid, k)
            uz, vz = τbx0 * (1 - ζ) / μ0, τby0 * (1 - ζ) / μ0
            expected = sqrt(horizontal + (uz^2 + vz^2) / 4)
            @test all(≈(expected, rtol = 1e-12), interior(mech.strainrate.effective)[:, :, k])
        end
        # The horizontal part alone would be depth-independent; confirm the layers
        # actually differ, i.e. the shear term is doing something rather than vanishing.
        eff = interior(mech.strainrate.effective)
        @test !(eff[4, 4, 1] ≈ eff[4, 4, grid.nz])
    end

    # f̄ = ∫₀¹ f dζ. Exact to roundoff for a uniform field on any layering (weights sum to
    # 1 by construction); to the midpoint rule's own O(Δζ²) once f varies with z, matching
    # the accuracy behaviour `viscosity_integrals!` already established for F_m.
    @testset "depthaverage!: exact for uniform, quadrature-order for varying" begin
        lay  = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 10))
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, lay)
        rt   = Runtime(grid)
        mech = MechanicState(grid)
        out  = Field(rt.arch, rt.grid2d, (Center(), Center(), Center()))

        μ0 = 1e5
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0)
        depthaverage!(out, mech.material.viscosity, rt)
        @test all(==(μ0), interior(out))    # roundoff-exact, not merely close

        # ζ (not ζ²) is linear, and the midpoint rule is exact for a linear integrand too
        # (the same fact the Stage 1 F₁ test relies on) — so this needs a genuinely
        # quadratic profile to see any error at all.
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0 * (1 + ζ^2))
        depthaverage!(out, mech.material.viscosity, rt)
        exact = μ0 * (1 + 1 / 3)
        @test all(≈(exact, rtol = 2e-3), interior(out))
        @test !all(≈(exact, rtol = 1e-10), interior(out))

        # Masked-out columns contribute nothing rather than a division artifact.
        grid_flat = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        rt_flat, mech_flat = Runtime(grid_flat), MechanicState(grid_flat)
        topo = TopographicState(grid_flat)
        setdata!(mech_flat.material.viscosity, μ0)
        setdata!(topo.mask.is_ice, false)
        out_flat = Field(rt_flat.arch, rt_flat.grid2d, (Center(), Center(), Center()))
        depthaverage!(out_flat, mech_flat.material.viscosity, rt_flat, IceMask(topo.mask.is_ice))
        @test all(==(0.0), interior(out_flat))
    end

    # β_eff = β/(1+βF₂) = 1/(1/β + F₂) — Eq. 19, and Eq. 20's frozen-bed limit as β → ∞.
    # The reciprocal form (decision 11) is what makes both limits ordinary IEEE arithmetic
    # rather than a branch on a "large β" threshold.
    @testset "beta_eff_diva!: Eq. 19, spatially varying, and both limits" begin
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0)
        rt   = Runtime(grid)
        mech = MechanicState(grid)
        be   = Field(rt.arch, rt.grid2d, (Center(), Center(), Center()))

        # A genuinely spatially-varying case, not just uniform constants — catches an
        # accidental broadcast-scalar shortcut that a uniform test could miss.
        fill_analytic!(mech.friction.beta, rt.grid2d, (x, y) -> 100.0 + 5.0x)
        fill_analytic!(mech.material.viscosity_integral_2, rt.grid2d, (x, y) -> 1e-3 + 2e-4y)
        beta_eff_diva!(be, mech.friction.beta, mech.material.viscosity_integral_2, rt)
        ref = analytic_like(be, rt.grid2d, (x, y) -> begin
            β, F2 = 100.0 + 5.0x, 1e-3 + 2e-4y
            β / (1 + β * F2)
        end)
        @test all(isapprox.(interior(be), ref; rtol = 1e-12))

        F2val = 1e-3
        setdata!(mech.material.viscosity_integral_2, F2val)
        for (β, expected) in ((1e3, 1e3 / (1 + 1e3 * F2val)),   # generic Eq. 19
                              (0.0, 0.0),                        # free slip
                              (Inf, 1 / F2val))                  # frozen bed, Eq. 20
            setdata!(mech.friction.beta, β)
            beta_eff_diva!(be, mech.friction.beta, mech.material.viscosity_integral_2, rt)
            @test all(≈(expected, rtol = 1e-12), interior(be))
        end

        # SSA is the F₂ = 0 limit: β_eff degrades to β unchanged, not to nonsense, if the
        # integrals were never populated.
        setdata!(mech.friction.beta, 4.2e13)
        setdata!(mech.material.viscosity_integral_2, 0.0)
        beta_eff_diva!(be, mech.friction.beta, mech.material.viscosity_integral_2, rt)
        @test all(≈(4.2e13, rtol = 1e-12), interior(be))
    end

    # DIVAViscosityContinuation: Glen's law (same closed form as GlenViscosityContinuation,
    # see its own test above this one) evaluated per layer from Eq. 13's effective strain
    # rate, log-space relaxed exactly as the depth-averaged version — then
    # viscosity_depthaveraged derived by depthaverage!, not computed independently.
    @testset "DIVAViscosityContinuation: matches the closed-form formula, per layer" begin
        lay  = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 6))
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, lay)
        rt   = Runtime(grid)
        mech = MechanicState(grid)

        # µ_old0 is Pa yr, like every other viscosity here. It used to be `3e14`, a Pa s
        # value, which made the shear term τ_b(1-ζ)/µ ~ 1e-12 against a horizontal strain
        # rate of ~4e-3 — so every layer came out with the same µ and the per-layer
        # ζ-dependence this testset exists to check was never actually exercised.
        A0, n, ε̇0, μ_old0 = 1e-16, 3.0, 1e-12, 3e5
        a, b, c, d = 2e-3, -1e-3, 5e-4, 3e-3
        τbx0, τby0 = 500.0, -200.0

        fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d, (x, y) -> a * x + b * y)
        fill_analytic!(mech.velocity.depthaverage_y, rt.grid2d, (x, y) -> c * x + d * y)
        fill_analytic!(mech.stress.base_x, rt.grid2d, (x, y) -> τbx0)
        fill_analytic!(mech.stress.base_y, rt.grid2d, (x, y) -> τby0)
        fill_analytic3d!(mech.material.rate_factor, rt.grid, (x, y, ζ) -> A0)
        depthaverage_velocitygradients!(mech.velocity, rt)

        horizontal = a^2 + d^2 + a * d + ((b + c) / 2)^2
        μ_raw_at(ζ) = begin
            uz, vz = τbx0 * (1 - ζ) / μ_old0, τby0 * (1 - ζ) / μ_old0
            eff = sqrt(horizontal + (uz^2 + vz^2) / 4)
            inv(2 * A0^(1 / n)) * sqrt(eff^2 + ε̇0^2)^((1 - n) / n)
        end

        setdata!(mech.material.viscosity, μ_old0)
        vc1 = DIVAViscosityContinuation(; n_glen = n, theta_mu = 1.0, strainrate_reg = ε̇0)
        update_viscosity!(mech, vc1, rt)
        for k in 1:grid.nz
            @test all(≈(μ_raw_at(zcenter(rt.grid, k)), rtol = 1e-10),
                     interior(mech.material.viscosity)[:, :, k])
        end

        # The shear term must actually move µ from layer to layer, or the loop above is
        # just checking one number nz times (which is what the old Pa s seed did).
        μ_z = interior(mech.material.viscosity)
        @test !(μ_z[4, 4, 1] ≈ μ_z[4, 4, grid.nz])

        setdata!(mech.material.viscosity, μ_old0)
        vc2 = DIVAViscosityContinuation(; n_glen = n, theta_mu = 0.2, strainrate_reg = ε̇0)
        update_viscosity!(mech, vc2, rt)
        for k in 1:grid.nz
            expected = exp(0.2 * log(μ_raw_at(zcenter(rt.grid, k))) + 0.8 * log(μ_old0))
            @test all(≈(expected, rtol = 1e-10), interior(mech.material.viscosity)[:, :, k])
        end

        # viscosity_depthaveraged is *derived*, so it must equal depthaverage! of the µ(z)
        # just written — not independently close to it.
        reference = Field(rt.arch, rt.grid2d, (Center(), Center(), Center()))
        depthaverage!(reference, mech.material.viscosity, rt)
        @test interior(mech.material.viscosity_depthaveraged) == interior(reference)
    end

    # diva_update!'s dependency order (µ → µ̄ → F₁/F₂ → β_eff) end-to-end: F₂ matches the
    # uniform-µ closed form H/(3µ) (the paper's own check, quoted at their Eq. 15/§3.1
    # line "F₂ = H/(3µ)"), and β_eff matches Eq. 19 built from that F₂.
    #
    # `NoViscosityContinuation` deliberately: this isolates the F₁/F₂/β_eff half of the
    # chain, given a fixed `µ(z)`, from the "how µ itself is derived" half already covered
    # by the `DIVAViscosityContinuation` testset above — with the real continuation active,
    # `update_viscosity!` would overwrite the prefilled `µ0` from the (zero) velocity
    # gradients this setup uses, and the test would silently stop checking what it says it
    # checks.
    @testset "diva_update!: full chain, uniform column" begin
        lay  = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 24))
        grid = StaggeredGrid(Float64, 8.0, 8.0, 1.0, 1.0, lay)
        rt   = Runtime(grid)
        mech = MechanicState(grid)
        H0, μ0, β0 = 1000.0, 1e5, 1e4

        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
        fill_analytic!(mech.friction.beta, rt.grid2d, (x, y) -> β0)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0)

        solver = PseudoTransientSolver(grid;
            viscosity_continuation = NoViscosityContinuation())
        diva_update!(mech, solver, rt)

        F2exact  = H0 / (3μ0)
        βeff_exact = β0 / (1 + β0 * F2exact)
        @test all(≈(F2exact, rtol = 1e-2), interior(mech.material.viscosity_integral_2))
        @test all(≈(βeff_exact, rtol = 1e-2), interior(mech.friction.beta_eff))
    end

    # Robinson et al. (2022) §3.1's own reference state, DIVA's simplest case: uniform
    # thickness/viscosity/slope/friction on a frozen (zero-velocity start) slab. The paper
    # gives F₂ = H/(3µ) explicitly (their line above Eq. 55) and Eq. 18 then gives the
    # depth-averaged velocity as sliding *plus* the vertical-shear correction —
    # ū = u_b(1+βF₂) = τ_d/β + τ_d·H/(3µ), strictly larger than SSA's τ_d/β alone. This is
    # the one place in the whole port where DIVA and SSA give genuinely different numbers.
    @testset "DIVA linearised slab ($_DIVA_PAPER, §3.1)" begin
        cst = Constants{Float64}()
        H0, μ0, β0, α = 1000.0, 1e5, 1e4, 1e-3
        dx = 5e3
        τd = cst.density_ice * cst.gravity * H0 * α
        ub_ssa  = τd / β0

        function setup(nz)
            lay  = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, nz))
            grid = StaggeredGrid(Float64, 10dx, 3dx, dx, dx, lay)
            rt   = Runtime(grid)
            mech = MechanicState(grid)
            fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
            fill_analytic!(mech.topography.surface, rt.grid2d, (x, y) -> H0 - α * x)
            fill_analytic!(mech.material.viscosity_depthaveraged, rt.grid2d, (x, y) -> μ0)
            fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0)
            fill_analytic!(mech.friction.beta, rt.grid2d, (x, y) -> β0)
            return grid, rt, mech
        end

        # Convergence in nz: the whole discrepancy from the exact ū is F₂'s O(Δζ²)
        # quadrature error (Stage 1), so it must shrink like the F_m tests already showed.
        errs = Float64[]
        for nz in (16, 32, 64)
            grid, rt, mech = setup(nz)
            F2exact = H0 / (3μ0)
            ū_exact = ub_ssa + τd * F2exact

            solver = PseudoTransientSolver(grid; maxiter = 500, abstol = 1e-8)
            diva_update!(mech, solver, rt)   # caller's responsibility under NoDIVUpdate
            res = pseudo_transient!(mech, cst, solver, rt, DIVAMomentumBalance())

            @test res.converged
            ū = interior(mech.velocity.depthaverage_x)
            @test all(≈(ū_exact, rtol = 5e-2), ū)   # loose: pins gross correctness first
            push!(errs, maximum(abs, ū .- ū_exact) / ū_exact)

            if nz == 32
                # The point of this whole exercise: DIVA sliding is not SSA sliding. β0/µ0
                # here are chosen so the shear correction is the *dominant* term (paper's
                # own regime, §3.1's "high aspect ratio" case) — a bug that silently fell
                # back to the SSA path would be off by an order of magnitude, not a
                # rounding error.
                @test all(>(10), ū ./ ub_ssa)
            end
        end
        @test all(>(1.7), convergence_rates(errs))   # ~2nd order, matching F₂'s quadrature
    end

    # The time-unit anchor. Every other slab test here prescribes µ̄ directly, which makes
    # them blind to the convention: µ and β scale together under a change of time unit,
    # τ_d is unaffected (Pa carries no time), and ū comes out right in *either* system.
    # This one starts from the rate factor instead, so `A`'s time unit and the velocity's
    # are forced to agree — a half-applied seconds/years conversion in the viscosity chain
    # shows up as a factor of ~3.16e7, not a rounding error.
    #
    # Pure vertical shear in a uniform slab (no membrane stress under periodic BCs) has a
    # closed form: τ_xz(z) = ρgα(H-z) and du/dz = 2A τ_xz^n integrate to
    #   ū = u_b + (2A/(n+2))·(ρgα)^n·H^(n+1),   u_b = τ_d/β.
    # β is deliberately stiff so the deformational term dominates sliding ~30×: the point
    # is to measure `A`, not to re-test τ_d/β, which the tests above already cover.
    #
    # The absolute number is what pins the convention for a human reader: A = 1e-16 is
    # Cuffey & Paterson's Pa⁻³ **yr⁻¹** value, and ~3 cm/yr is a sane interior-ice-sheet
    # speed for a 1 km slab on a 0.1% slope. Read as Pa⁻³s⁻¹ it would be ~9e5 m/yr.
    @testset "Glen-law slab from the rate factor pins the time unit" begin
        cst = Constants{Float64}()
        H0, β0, α = 1000.0, 1e7, 1e-3
        A0, n = 1e-16, 3.0
        dx = 5e3

        ρg = cst.density_ice * cst.gravity
        τd = ρg * H0 * α
        ū_exact = τd / β0 + (2A0 / (n + 2)) * (ρg * α)^n * H0^(n + 1)

        # `A` sets the viscosity scale, so seed µ with the value it implies at the bed
        # rather than an unrelated constant — a wildly wrong start just costs iterations.
        μ_seed = 1 / (2 * A0 * τd^(n - 1))

        function setup(nz)
            lay  = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, nz))
            grid = StaggeredGrid(Float64, 10dx, 3dx, dx, dx, lay)
            rt   = Runtime(grid)
            mech = MechanicState(grid)
            fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
            fill_analytic!(mech.topography.surface, rt.grid2d, (x, y) -> H0 - α * x)
            fill_analytic!(mech.material.viscosity_depthaveraged, rt.grid2d, (x, y) -> μ_seed)
            fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ_seed)
            fill_analytic3d!(mech.material.rate_factor, rt.grid, (x, y, ζ) -> A0)
            fill_analytic!(mech.friction.beta, rt.grid2d, (x, y) -> β0)
            return grid, rt, mech
        end

        for nz in (32, 64)
            grid, rt, mech = setup(nz)
            # ε̇₀ must sit well below the slab's own strain rate (~1.4e-4 yr⁻¹ at the bed)
            # or the floor, not Glen's law, would set the viscosity.
            vc = DIVAViscosityContinuation(; n_glen = n, theta_mu = 1.0,
                                             strainrate_reg = 1e-10)
            solver = PseudoTransientSolver(grid; maxiter = 20_000, abstol = 1e-12,
                                           viscosity_continuation = vc,
                                           div_update = PeriodicDIVUpdate(1))
            diva_update!(mech, solver, rt)
            res = pseudo_transient!(mech, cst, solver, rt, DIVAMomentumBalance())

            @test res.converged
            # Domain centre only. The closed form is for an unbounded slab, and `surface`
            # here is a linear ramp, which cannot be periodic — so the wrap leaves the end
            # x-columns with a wrong surface gradient. With a *fixed* viscosity (every
            # other slab test above) that error cannot feed back and the field stays
            # uniform right to the edge; under a nonlinear viscosity the bad edge strain
            # rate sets a bad edge viscosity, and the error diffuses inward over ~3 cells
            # (halved at the edge itself, ~7% one cell in, under 1% by the centre). That
            # boundary layer belongs to the ramp, not to the rheology being measured.
            ū = interior(mech.velocity.depthaverage_x)[4:8, :, :]
            @test all(≈(ū_exact, rtol = 3e-2), ū)
        end
    end

    # The rate factors must agree with each other on what a second is. `A(T) = A₀exp(-Q/RT)`
    # has Q, R and T all time-free, so A₀ alone carries the unit — which makes the
    # temperature-dependent laws directly comparable against `PrescribedRateFactor`'s
    # constant with no free scale to hide a conversion error. Before the seconds→years
    # migration these disagreed by ~2e8, i.e. exactly `SECONDS_PER_YEAR`, and swapping one
    # for the other through the shared `AbstractRateFactor` interface froze the ice solid.
    @testset "rate factors agree on the time unit" begin
        A_ref = PrescribedRateFactor().A          # 1e-16 Pa⁻³ yr⁻¹, Cuffey & Paterson

        # Near the melting point the piecewise-Arrhenius law should land within a small
        # factor of the constant everyone quotes. A unit slip is eight orders, not two.
        for rf in (ArrheniusRateFactor(), LliboutryDuvalRateFactor())
            A = rate_factor(272.15, rf)
            @test 0.1 < A / A_ref < 10
        end

        # Hooke gets a looser band: it runs ~250x above Arrhenius at *every* temperature,
        # not just near melting, so the offset is in `A_0` rather than the
        # proximity-to-melting term. That is a pre-existing calibration question — the
        # ratio was identical when both prefactors were still in seconds — and this test
        # is only here to catch unit slips, which are five orders larger than the gap.
        @test 1 < rate_factor(272.15, HookeRateFactor()) / A_ref < 1e4

        # Colder ice is stiffer, and monotonically so.
        A_cold = rate_factor(253.15, ArrheniusRateFactor())
        A_warm = rate_factor(272.15, ArrheniusRateFactor())
        @test A_cold < A_warm

        # `time_unit = :second` must undo exactly the conversion baked into the defaults.
        @test ArrheniusRateFactor(:second; A_0_p1 = 3.985e-13, A_0_p2 = 1.916e3) ==
              ArrheniusRateFactor()
        @test HookeRateFactor(:second; A_0 = 9.302e-7) == HookeRateFactor()
        @test PrescribedRateFactor(:second; A = 3.2e-24).A ≈ A_ref rtol = 0.02
        @test PrescribedRateFactor(:year; A = 1e-16) == PrescribedRateFactor()

        # Time-free fields must pass through untouched, and partial overrides must leave
        # the other prefactors at their (already internal) defaults.
        @test ArrheniusRateFactor(:second; Q_a_p1 = 60e3).Q_a_p1 == 60e3
        @test ArrheniusRateFactor(:second; A_0_p1 = 3.985e-13) == ArrheniusRateFactor()
        @test_throws ArgumentError ArrheniusRateFactor(:fortnight)
    end

    # PeriodicDIVUpdate must actually reach `pseudo_dt!`, not just `diva_update!` — decision
    # 8's whole point is that a grown β_eff under a stale Δτ bound diverges.
    #
    # Comparing `dtau_x` before vs. after one `pseudo_transient!` call is the wrong test:
    # `_tuning_init!` already calls `pseudo_dt!` once unconditionally before the loop even
    # starts, for *every* strategy, so `dtau_x` is never literally "untouched". The real
    # invariant is whether `dtau_x` is a *function of how many iterations ran*: under
    # `NoDIVUpdate` nothing inside the loop can change the Gershgorin bound's inputs (no
    # continuation touches `β_eff`/`µ̄` after that one initial call), so `maxiter = 1` and
    # `maxiter = 5` must land on bit-identical `dtau_x`. Under `PeriodicDIVUpdate(1)` with
    # `DIVAViscosityContinuation` active, `µ(z)` genuinely evolves with the velocity
    # iterate each PT step, so the two must differ.
    @testset "dtau_x depends on iteration count under PeriodicDIVUpdate, not NoDIVUpdate" begin
        lay = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 16))
        dx = 5e3
        H0, μ0, β0, α = 1000.0, 1e5, 1e4, 1e-3
        cst = Constants{Float64}()

        function setup()
            grid = StaggeredGrid(Float64, 10dx, 3dx, dx, dx, lay)
            rt   = Runtime(grid)
            mech = MechanicState(grid)
            fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
            fill_analytic!(mech.topography.surface, rt.grid2d, (x, y) -> H0 - α * x)
            fill_analytic!(mech.material.viscosity_depthaveraged, rt.grid2d, (x, y) -> μ0)
            fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0)
            fill_analytic!(mech.friction.beta, rt.grid2d, (x, y) -> β0)
            fill_analytic3d!(mech.material.rate_factor, rt.grid, (x, y, ζ) -> 1e-16)
            return grid, rt, mech
        end

        function run(maxiter, div_update, vc)
            grid, rt, mech = setup()
            solver = PseudoTransientSolver(grid; maxiter, div_update,
                                           viscosity_continuation = vc)
            diva_update!(mech, solver, rt)
            pseudo_transient!(mech, cst, solver, rt, DIVAMomentumBalance())
            return copy(interior(solver.dtau_x))
        end

        vc = NoViscosityContinuation()
        dtau_no_1 = run(1, NoDIVUpdate(), vc)
        dtau_no_5 = run(5, NoDIVUpdate(), vc)
        @test dtau_no_1 == dtau_no_5

        vc_diva = DIVAViscosityContinuation(; strainrate_reg = 1e-12)
        dtau_periodic_1 = run(1, PeriodicDIVUpdate(1), vc_diva)
        dtau_periodic_5 = run(5, PeriodicDIVUpdate(1), vc_diva)
        @test dtau_periodic_1 != dtau_periodic_5

        @test_throws ArgumentError PeriodicDIVUpdate(0)
    end
end

using Pagos
using Test

include("../test_helpers/chmy.jl")

# The Chmy-native DIVA viscosity integrals `F_m = ∫_b^s (1/µ)((s-z)/H)^m dz`
# (Robinson et al. 2022, Eq. 15). Validated against closed-form σ-integrals, per the hybrid
# migration strategy — never against the collocated `aggregated_viscosity_integral!`, which
# takes its σ levels from a caller-supplied vector rather than from the grid.
#
# Every case below uses a *stretched* (quadratic) layering, not a uniform one: the whole
# point of reading `zcenter`/`Δz` off the grid instead of taking a `sigma` argument is that
# the midpoint rule stays exact-by-construction when the layers are not evenly spaced, and
# a uniform layering would let a wrong-but-symmetric quadrature pass.

# F_m under a viscosity that varies with σ as µ(σ) = µ0 / (1 + σ):
#   F_1 = (H/µ0) ∫₀¹ (1-σ)(1+σ) dσ = 2H  / (3µ0)
#   F_2 = (H/µ0) ∫₀¹ (1-σ)²(1+σ) dσ = 5H / (12µ0)
const _F1_VARYING = 2 / 3
const _F2_VARYING = 5 / 12

function setup_integrals(T, nz; lx = 4.0, ly = 4.0, dx = 1.0, dy = 1.0)
    layering = CorrectedVerticalLayering(T, QuadraticSigmaTransform(T, nz))
    grid = StaggeredGrid(T, T(lx), T(ly), T(dx), T(dy), layering)
    rt   = Runtime(grid)
    mech = MechanicState(grid)
    F1 = Field(rt.arch, rt.grid2d, (Center(), Center(), Center()))
    F2 = Field(rt.arch, rt.grid2d, (Center(), Center(), Center()))
    return grid, rt, mech, F1, F2
end

@testset "DIVA viscosity integrals (C-grid staggered)" begin

    @testset "uniform µ: F₁ exact on a stretched layering, F₂ to quadrature order" begin
        H0, μ0 = 1000.0, 1e5
        grid, rt, mech, F1, F2 = setup_integrals(Float64, 12)

        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0)

        viscosity_integrals!(F1, F2, mech, rt)

        # The midpoint rule integrates a linear integrand exactly, and (1-σ)/µ is linear
        # under uniform µ — so this holds to roundoff no matter how stretched the layers
        # are. It is the assertion that pins "σ comes from the grid": a quadrature that
        # used interfaces, or midpoints of the wrong axis, would miss here.
        @test all(≈(H0 / (2μ0), rtol = 1e-14), interior(F1))

        # (1-σ)²/µ is quadratic, so F₂ carries the midpoint rule's O(Δζ²) error and is
        # only close.
        @test all(≈(H0 / (3μ0), rtol = 1e-2), interior(F2))
        @test !all(≈(H0 / (3μ0), rtol = 1e-14), interior(F2))
    end

    @testset "depth-varying µ: both integrals second-order convergent" begin
        H0, μ0 = 1000.0, 1e5
        err1, err2 = Float64[], Float64[]

        for nz in (8, 16, 32)
            grid, rt, mech, F1, F2 = setup_integrals(Float64, nz)
            fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
            fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0 / (1 + ζ))

            viscosity_integrals!(F1, F2, mech, rt)

            push!(err1, maximum(abs, interior(F1) .- _F1_VARYING * H0 / μ0))
            push!(err2, maximum(abs, interior(F2) .- _F2_VARYING * H0 / μ0))
        end

        # Halving the layer thickness quarters the error: the midpoint rule is second
        # order, and nothing in the σ → z scaling degrades it.
        @test all(>(1.9), convergence_rates(err1))
        @test all(>(1.9), convergence_rates(err2))
        @test err1[end] / (_F1_VARYING * H0 / μ0) < 1e-3
        @test err2[end] / (_F2_VARYING * H0 / μ0) < 1e-3
    end

    @testset "F_m scales linearly with H and inversely with µ" begin
        μ0 = 1e5
        grid, rt, mech, F1, F2 = setup_integrals(Float64, 10)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0 / (1 + ζ))

        # H varying in x, µ independent of it ⟹ F_m(x) / H(x) is constant.
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> 500.0 + 100.0 * x)
        viscosity_integrals!(F1, F2, mech, rt)
        Hi = analytic_like(mech.topography.thickness, rt.grid2d, (x, y) -> 500.0 + 100.0 * x)
        ratio1, ratio2 = interior(F1) ./ Hi, interior(F2) ./ Hi
        @test all(≈(first(ratio1), rtol = 1e-13), ratio1)
        @test all(≈(first(ratio2), rtol = 1e-13), ratio2)

        # Doubling µ everywhere halves both integrals.
        F1b = Field(rt.arch, rt.grid2d, (Center(), Center(), Center()))
        F2b = Field(rt.arch, rt.grid2d, (Center(), Center(), Center()))
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> 2μ0 / (1 + ζ))
        viscosity_integrals!(F1b, F2b, mech, rt)
        @test all(≈(0.5, rtol = 1e-13), interior(F1b) ./ interior(F1))
        @test all(≈(0.5, rtol = 1e-13), interior(F2b) ./ interior(F2))
    end

    @testset "ice-free column: H = 0 gives 0, not a division by zero" begin
        μ0 = 1e5
        grid, rt, mech, F1, F2 = setup_integrals(Float64, 10)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> x < 0 ? 0.0 : 800.0)

        viscosity_integrals!(F1, F2, mech, rt)

        @test !any(isnan, interior(F1))
        @test !any(isinf, interior(F1))
        for i in axes(interior(F1), 1)
            x, _, _ = coord(rt.grid2d, location(F1), i, 1, 1)
            expected = x < 0 ? 0.0 : 800.0 / (2μ0)
            @test all(≈(expected, rtol = 1e-13), interior(F1)[i, :, 1])
        end
    end

    # Containment: the integrand is 1/µ, so an unmasked ice-free cell carrying a zero
    # placeholder viscosity produces Inf — mirroring the `hlerp`-gives-NaN pattern the
    # membrane-stress tests already establish. Note H > 0 here, so the H guard cannot help;
    # only the mask can.
    @testset "µ = 0 gives Inf off-ice, contained by IceMask" begin
        grid, rt, mech, F1, F2 = setup_integrals(Float64, 10)
        topo = TopographicState(grid)

        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> 800.0)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> x < 0 ? 0.0 : 1e5)

        viscosity_integrals!(F1, F2, mech, rt)
        @test any(isinf, interior(F1))
        @test any(isinf, interior(F2))

        setdata!(topo.mask.is_ice, true)
        for i in axes(interior(topo.mask.is_ice), 1)
            x, _, _ = coord(rt.grid2d, location(topo.mask.is_ice), i, 1, 1)
            x < 0 && (interior(topo.mask.is_ice)[i, :, 1] .= false)
        end

        viscosity_integrals!(F1, F2, mech, rt, IceMask(topo.mask.is_ice))
        @test !any(isinf, interior(F1))
        @test !any(isinf, interior(F2))
        @test !any(isnan, interior(F1))
    end

    # The documented `nz == 1` limitation, pinned so it cannot regress silently into
    # "looks fine": one layer puts the quadrature point at σ = 1/2, which is exact for the
    # linear F₁ integrand and 25% low for the quadratic F₂ one.
    @testset "nz == 1 resolves F₁ but not F₂" begin
        H0, μ0 = 1000.0, 1e5
        grid = StaggeredGrid(Float64, 4.0, 4.0, 1.0, 1.0)
        rt   = Runtime(grid)
        mech = MechanicState(grid)
        @test rt.grid === rt.grid2d

        F1 = Field(rt.arch, rt.grid2d, (Center(), Center(), Center()))
        F2 = Field(rt.arch, rt.grid2d, (Center(), Center(), Center()))
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0)

        viscosity_integrals!(F1, F2, mech, rt)

        @test all(≈(H0 / (2μ0), rtol = 1e-14), interior(F1))
        @test all(≈(H0 / (4μ0), rtol = 1e-14), interior(F2))   # not H0/(3μ0)
    end

    @testset "Float32 stays Float32" begin
        H0, μ0 = 1000f0, 1f5
        grid, rt, mech, F1, F2 = setup_integrals(Float32, 12)

        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)
        fill_analytic3d!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0)

        viscosity_integrals!(F1, F2, mech, rt)

        @test eltype(F1) === Float32
        @test eltype(F2) === Float32
        @test all(≈(H0 / (2μ0), rtol = 1f-6), interior(F1))
    end
end

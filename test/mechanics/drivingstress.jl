using Pagos
using Test

include("../test_helpers/chmy.jl")

# The Chmy-native, C-grid staggered driving stress: `τ = ρ g H ∂s/∂x` with the gradient
# taken natively `aa → acx`/`acy` and `H` lerped onto the same face. Validated against
# analytic solutions, not against the collocated method it replaces (which discretizes a
# different thing — a cell-centred gradient — see `roadmaps/chmy.md`, Phase 3).

@testset "driving stress (C-grid staggered)" begin
    lx, ly, dx, dy = 8.0, 8.0, 1.0, 1.0
    cst = Constants{Float64}()
    ρg  = cst.density_ice * cst.gravity

    function setup(T = Float64)
        grid = StaggeredGrid(T, lx, ly, dx, dy)
        return grid, Runtime(grid), MechanicState(grid), TopographicState(grid)
    end

    # The whole point of the port: the driving stress must land where the velocity it
    # forces lives, with no interpolation of the gradient.
    @testset "the driving stress lands on the velocity faces" begin
        grid, _, mech, _ = setup()
        @test location(mech.stress.driving_x) === location(mech.velocity.x_bar)
        @test location(mech.stress.driving_y) === location(mech.velocity.y_bar)
        @test size(interior(mech.stress.driving_x)) == (grid.nx + 1, grid.ny, 1)
        @test size(interior(mech.stress.driving_y)) == (grid.nx, grid.ny + 1, 1)
    end

    # Linear surface, uniform thickness: `lerp` of a constant is exact and the two-point
    # gradient of a linear function is exact, so the result is exact everywhere.
    @testset "exact: linear surface, uniform H" begin
        a, b, H0 = 1000.0, 0.02, 800.0
        _, rt, mech, _ = setup()

        fill_analytic!(mech.topography.surface, rt.grid2d, (x, y) -> a + b * x)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)

        drivingstress!(mech, cst, rt)

        @test all(interior(mech.stress.driving_x) .≈ ρg * H0 * b)
        @test all(interior(mech.stress.driving_y) .≈ 0.0)
    end

    # Same in y, to catch an axis mix-up: a `Dim(1)` left in the y branch, or `NODE_ACX`
    # used for the y face, would pass the x-only test above.
    @testset "exact: linear surface in y, uniform H" begin
        a, b, H0 = 500.0, -0.01, 1200.0
        _, rt, mech, _ = setup()

        fill_analytic!(mech.topography.surface, rt.grid2d, (x, y) -> a + b * y)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)

        drivingstress!(mech, cst, rt)

        @test all(interior(mech.stress.driving_x) .≈ 0.0)
        @test all(interior(mech.stress.driving_y) .≈ ρg * H0 * b)
    end

    # A stronger exact case, because both factors vary and both must be staggered
    # correctly. A two-point difference of a *quadratic* is exact at the face it spans
    # (`(x_i² - x_{i-1}²)/Δx == x_i + x_{i-1} == 2·x_face`), and `lerp` of a linear H is
    # exact — so the product is exact at every face, with no truncation error to hide a
    # half-cell shift. Getting either stagger wrong moves the answer off by O(Δx).
    @testset "exact: quadratic surface, linear H" begin
        c, a, b = 3.0e-4, 900.0, 5.0
        _, rt, mech, _ = setup()

        fill_analytic!(mech.topography.surface, rt.grid2d, (x, y) -> c * x^2)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> a + b * x)

        drivingstress!(mech, cst, rt)

        expected = analytic_like(mech.stress.driving_x, rt.grid2d,
                                (x, y) -> ρg * (a + b * x) * 2c * x)
        @test interior(mech.stress.driving_x) ≈ expected
    end

    # The textbook uniform-slab result every ice-sheet model is checked against:
    # |τ_d| = ρ g H tan α on a plane of slope α.
    @testset "uniform slab: |τ_d| = ρ g H tanα" begin
        α, H0 = deg2rad(0.5), 1000.0
        _, rt, mech, _ = setup()

        fill_analytic!(mech.topography.surface, rt.grid2d, (x, y) -> -tan(α) * x)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)

        drivingstress!(mech, cst, rt)
        @test all(abs.(interior(mech.stress.driving_x)) .≈ ρg * H0 * tan(α))
    end

    # The sign the pseudo-transient residual depends on: `dotvel!` *subtracts* this field,
    # so what is stored is `+ρgH∇s`, not the physical `-ρgH∇s`. Flipping it would silently
    # reverse the flow direction of every momentum balance.
    @testset "sign convention: stores +ρgH∇s" begin
        _, rt, mech, _ = setup()
        fill_analytic!(mech.topography.surface, rt.grid2d, (x, y) -> 0.02x)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> 800.0)

        drivingstress!(mech, cst, rt)
        # surface rising with x ⟹ stored component positive (ice flows toward -x)
        @test all(interior(mech.stress.driving_x) .> 0)
    end

    # Second-order convergence where the discretization is *not* exact.
    @testset "manufactured solution: 2nd-order convergence" begin
        k = 2π / lx
        sfun(x, y) = 100 * sin(k * x) * cos(k * y)
        Hfun(x, y) = 1000 + 200 * cos(k * x)
        τxfun(x, y) = ρg * Hfun(x, y) * 100 * k * cos(k * x) * cos(k * y)

        errors = map((16, 32, 64, 128)) do n
            h    = lx / n
            grid = StaggeredGrid(Float64, lx, ly, h, h)
            rt   = Runtime(grid)
            mech = MechanicState(grid)

            fill_analytic!(mech.topography.surface, rt.grid2d, sfun)
            fill_analytic!(mech.topography.thickness, rt.grid2d, Hfun)
            drivingstress!(mech, cst, rt)

            exact = analytic_like(mech.stress.driving_x, rt.grid2d, τxfun)
            maximum(abs, interior(mech.stress.driving_x) .- exact) / (ρg * 1000 * 100 * k)
        end

        @test issorted(errors; rev = true)
        @test all(r -> 1.9 < r < 2.1, convergence_rates(errors))
    end

    # `surface_gradient!` and `drivingstress!` compute the gradient by separate code paths
    # (the latter inline, to stay a single sweep), so they must agree — and the stored
    # gradients must be at the faces, not the centres.
    @testset "surface_gradient! agrees with the driving stress' inline gradient" begin
        _, rt, mech, topo = setup()
        sfun(x, y) = 0.01x^2 - 0.02y + 3
        H0 = 700.0

        for f in (mech.topography.surface, topo.elevation.surface)
            fill_analytic!(f, rt.grid2d, sfun)
        end
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> H0)

        surface_gradient!(topo, rt)
        drivingstress!(mech, cst, rt)

        @test location(topo.elevation.surface_dx) === location(mech.stress.driving_x)
        @test interior(mech.stress.driving_x) ≈
              ρg * H0 .* interior(topo.elevation.surface_dx)
        @test interior(mech.stress.driving_y) ≈
              ρg * H0 .* interior(topo.elevation.surface_dy)

        # ...and the gradient itself is right: ∂/∂x of 0.01x² is 0.02x, exact at the face.
        @test interior(topo.elevation.surface_dx) ≈
              analytic_like(topo.elevation.surface_dx, rt.grid2d, (x, y) -> 0.02x)
        @test all(interior(topo.elevation.surface_dy) .≈ -0.02)
    end

    # `ρ_ice * g` is a Float64 product of Float32 constants unless it is converted, which
    # would silently promote the whole kernel back to Float64.
    @testset "Float32 stays Float32" begin
        cst32 = Constants{Float32}()
        _, rt, mech, _ = setup(Float32)

        fill_analytic!(mech.topography.surface, rt.grid2d, (x, y) -> 0.02f0 * x)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> 800.0f0)

        drivingstress!(mech, cst32, rt)

        @test eltype(mech.stress.driving_x) === Float32
        @test eltype(interior(mech.stress.driving_x)) === Float32
        @test all(interior(mech.stress.driving_x) .≈
                  cst32.density_ice * cst32.gravity * 800.0f0 * 0.02f0)
    end

    # The one place where applying a collocated kernel to a Field-based state would be
    # silently wrong rather than an error: `deviatoric_stress!`'s fused flat index assumes
    # all tensor components share a shape, and on the C-grid they do not (four node
    # classes, three lengths, and the shortest one drives `ndrange`, so nothing throws).
    @testset "deviatoric_stress! rejects a Field-based state" begin
        layering = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 6))
        grid = StaggeredGrid(Float64, lx, ly, dx, dy, layering)
        mech, mat = MechanicState(grid), MaterialState(grid)

        # the premise of the guard: the components genuinely differ in shape
        @test size(mech.stress.xx) != size(mech.stress.xy)
        @test length(mech.stress.xx) < length(mech.stress.xy)

        @test_throws ErrorException deviatoric_stress!(mech, mat)

        # ...while the plain-array state still works, unchanged.
        rmech = MechanicState(RegularGrid(Float64, lx, ly, dx, dy))
        rmat  = MaterialState(RegularGrid(Float64, lx, ly, dx, dy))
        @test deviatoric_stress!(rmech, rmat) === nothing
    end

    # The collocated methods are untouched and still dispatch on their own signature
    # (8 positional args ending in dx, dy) rather than being shadowed by the Runtime one.
    @testset "the collocated method still exists" begin
        grid = RegularGrid(Float64, lx, ly, dx, dy)
        (; nx, ny) = grid
        τx, τy = zeros(nx, ny), zeros(nx, ny)
        s = [0.02 * x for x in grid.x, _ in grid.y]
        H = fill(800.0, nx, ny)

        drivingstress!(τx, τy, s, H, cst.density_ice, cst.gravity, dx, dy)
        @test τx[3, 3] ≈ ρg * 800.0 * 0.02
    end
end

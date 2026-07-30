using Pagos
using Test

include("../test_helpers/chmy.jl")

# Mass continuity, `∂H/∂t = -∇·q + ṁ`, validated against analytic/manufactured solutions
# rather than against the pre-migration code path (see `roadmaps/chmy.md`, Phase 3).

# Total mass rate over the interior, and the net inflow across the interior's outer faces.
# For a flux-form scheme these must agree with the total mass balance to machine
# precision — the interior flux differences telescope away exactly.
total_rate(dHdt, dx, dy) = sum(interior(dHdt)) * dx * dy

function boundary_influx(q_x, q_y, dx, dy)
    qx = interior(q_x)
    qy = interior(q_y)
    # inflow at the west/south faces of the interior block minus outflow at east/north
    return (sum(qx[1, :, :]) - sum(qx[end, :, :])) * dy +
           (sum(qy[:, 1, :]) - sum(qy[:, end, :])) * dx
end

@testset "advection (mass continuity)" begin
    lx, ly, dx, dy = 8.0, 8.0, 1.0, 1.0

    function setup(T = Float64; nz = nothing)
        grid = isnothing(nz) ? StaggeredGrid(T, lx, ly, dx, dy) :
               StaggeredGrid(T, lx, ly, dx, dy,
                             CorrectedVerticalLayering(T, QuadraticSigmaTransform(T, nz)))
        rt   = Runtime(grid)
        topo = TopographicState(grid)
        mech = MechanicState(grid)
        return grid, rt, topo, mech
    end

    @testset "state layout: fluxes are on the velocity faces" begin
        grid, _, _, mech = setup()

        # The Phase-2 deferred item: `flux` was a single `aa` scalar, which the divergence
        # would have had to re-stagger. It is now a face pair at the same node classes as
        # the depth-averaged velocity it is built from.
        @test location(mech.flux.x) === location(mech.velocity.x_bar)
        @test location(mech.flux.y) === location(mech.velocity.y_bar)
        @test size(interior(mech.flux.x)) == (grid.nx + 1, grid.ny, 1)
        @test size(interior(mech.flux.y)) == (grid.nx, grid.ny + 1, 1)
        @test size(interior(mech.flux.grline)) == (grid.nx, grid.ny, 1)

        # ...and the divergence of that pair lands exactly on the thickness field.
        @test location(mech.flux.grline) === location(mech.topography.thickness)
    end

    # `H` linear in x and `u` uniform. The centred average of a linear function is exact,
    # so `q_x = u·H(x_face)`; upwind instead reads the cell west of the face (`u > 0`), so
    # its face value is offset by half a cell — a genuinely different flux. The *tendency*
    # is nonetheless exact for both, because a uniform half-cell offset of a linear H
    # cancels in the difference: `∂q_x/∂x = u·b` either way.
    @testset "exact: linear H, uniform ū" begin
        a, b, U, V = 100.0, 3.0, 2.0, 0.0

        for scheme in (CenteredAdvection(), UpwindAdvection())
            grid, rt, topo, mech = setup()
            H, u, v = topo.thickness.ice, mech.velocity.x_bar, mech.velocity.y_bar

            fill_analytic!(H, rt.grid2d, (x, y) -> a + b * x)
            fill_analytic!(u, rt.grid2d, (x, y) -> U)
            fill_analytic!(v, rt.grid2d, (x, y) -> V)

            advect!(topo, mech, scheme, rt)

            xv     = collect(vertices(rt.grid2d, Dim(1)))
            xface  = scheme isa UpwindAdvection ? xv .- dx / 2 : xv
            @test interior(mech.flux.x)[:, 1, 1] ≈ U .* (a .+ b .* xface)
            @test all(interior(mech.flux.y) .≈ 0.0)
            @test all(interior(topo.thickness.ice_dt) .≈ -U * b)
        end
    end

    # Same in y, and with a mass-balance source, to catch an axis mix-up (a `Dim(1)` left
    # in the y branch would pass the x-only test above).
    @testset "exact: linear H in y, uniform v̄, with mass balance" begin
        a, b, V, mb = 50.0, -1.5, 4.0, 7.0
        grid, rt, topo, mech = setup()
        H, u, v = topo.thickness.ice, mech.velocity.x_bar, mech.velocity.y_bar

        fill_analytic!(H, rt.grid2d, (x, y) -> a + b * y)
        fill_analytic!(u, rt.grid2d, (x, y) -> 0.0)
        fill_analytic!(v, rt.grid2d, (x, y) -> V)
        fill_analytic!(topo.massbalance.net, rt.grid2d, (x, y) -> mb)

        advect!(topo, mech, CenteredAdvection(), rt)

        yv = collect(vertices(rt.grid2d, Dim(2)))
        @test interior(mech.flux.y)[1, :, 1] ≈ V .* (a .+ b .* yv)
        @test all(interior(topo.thickness.ice_dt) .≈ mb - V * b)
    end

    # Divergent flow: both components vary, so the two flux differences must both be right
    # *and* be added. q = H ū with H = const and ū = (αx, βy) gives ∇·q = H(α + β).
    @testset "exact: uniform H, linearly divergent ū" begin
        H0, α, β = 1000.0, 0.25, -0.1
        grid, rt, topo, mech = setup()

        fill_analytic!(topo.thickness.ice, rt.grid2d, (x, y) -> H0)
        fill_analytic!(mech.velocity.x_bar, rt.grid2d, (x, y) -> α * x)
        fill_analytic!(mech.velocity.y_bar, rt.grid2d, (x, y) -> β * y)

        advect!(topo, mech, CenteredAdvection(), rt)
        @test all(interior(topo.thickness.ice_dt) .≈ -H0 * (α + β))
    end

    @testset "upwind picks the upstream cell" begin
        grid, rt, topo, mech = setup()
        H, u, v = topo.thickness.ice, mech.velocity.x_bar, mech.velocity.y_bar

        # A thickness that is not symmetric about any face, so left ≠ right ≠ average.
        fill_analytic!(H, rt.grid2d, (x, y) -> exp(x / lx))
        fill_analytic!(v, rt.grid2d, (x, y) -> 0.0)

        Hi = interior(H)[:, 1, 1]
        i  = 3      # an interior x-face, with cells i-1 and i on either side

        fill_analytic!(u, rt.grid2d, (x, y) -> 1.0)     # eastward: upstream is cell i-1
        mass_flux!(mech.flux.x, mech.flux.y, H, u, v, UpwindAdvection(), rt)
        @test interior(mech.flux.x)[i, 1, 1] ≈ Hi[i - 1]

        fill_analytic!(u, rt.grid2d, (x, y) -> -1.0)    # westward: upstream is cell i
        mass_flux!(mech.flux.x, mech.flux.y, H, u, v, UpwindAdvection(), rt)
        @test interior(mech.flux.x)[i, 1, 1] ≈ -Hi[i]

        fill_analytic!(u, rt.grid2d, (x, y) -> 1.0)
        mass_flux!(mech.flux.x, mech.flux.y, H, u, v, CenteredAdvection(), rt)
        @test interior(mech.flux.x)[i, 1, 1] ≈ (Hi[i - 1] + Hi[i]) / 2
        @test !(interior(mech.flux.x)[i, 1, 1] ≈ Hi[i - 1])     # the schemes really differ
    end

    # The property the C-grid flux form exists for: total mass change = net boundary influx
    # + total mass balance, to *machine* precision (interior flux differences telescope),
    # for an arbitrary rough H and ū where the scheme is nowhere exact.
    @testset "discrete conservation (exact, not to truncation order)" begin
        mb = 0.3
        for scheme in (CenteredAdvection(), UpwindAdvection())
            grid, rt, topo, mech = setup()

            fill_analytic!(topo.thickness.ice, rt.grid2d,
                           (x, y) -> 500 + 100 * sin(3x) * cos(2y))
            fill_analytic!(mech.velocity.x_bar, rt.grid2d, (x, y) -> sin(x) + 0.5cos(y))
            fill_analytic!(mech.velocity.y_bar, rt.grid2d, (x, y) -> cos(2x) - 0.3sin(y))
            fill_analytic!(topo.massbalance.net, rt.grid2d, (x, y) -> mb)

            advect!(topo, mech, scheme, rt)

            expected = boundary_influx(mech.flux.x, mech.flux.y, dx, dy) +
                       mb * grid.nx * grid.ny * dx * dy
            @test total_rate(topo.thickness.ice_dt, dx, dy) ≈ expected rtol = 1e-13
        end
    end

    # Second-order convergence of the centred scheme on a smooth manufactured solution:
    # H = 1 + ½sin(kx), ū = (U, 0) ⟹ ∂H/∂t = -U·(k/2)cos(kx), analytically.
    @testset "manufactured solution: centred flux converges at 2nd order" begin
        U, k = 1.7, 2π / lx
        Hfun(x, y)    = 1 + 0.5 * sin(k * x)
        dHdtfun(x, y) = -U * 0.5 * k * cos(k * x)

        errors = map((16, 32, 64, 128)) do n
            h    = lx / n
            grid = StaggeredGrid(Float64, lx, ly, h, h)
            rt   = Runtime(grid)
            topo, mech = TopographicState(grid), MechanicState(grid)

            fill_analytic!(topo.thickness.ice, rt.grid2d, Hfun)
            fill_analytic!(mech.velocity.x_bar, rt.grid2d, (x, y) -> U)
            fill_analytic!(mech.velocity.y_bar, rt.grid2d, (x, y) -> 0.0)

            advect!(topo, mech, CenteredAdvection(), rt)

            exact = analytic_like(topo.thickness.ice_dt, rt.grid2d, dHdtfun)
            maximum(abs, interior(topo.thickness.ice_dt) .- exact)
        end

        @test issorted(errors; rev = true)
        rates = convergence_rates(errors)
        @test all(r -> 1.9 < r < 2.1, rates)
    end

    # First-order upwind on the same problem: converges, but only at order 1. Pins that
    # the two schemes are genuinely different discretizations and not the same code path.
    @testset "manufactured solution: upwind converges at 1st order" begin
        U, k = 1.7, 2π / lx
        Hfun(x, y)    = 1 + 0.5 * sin(k * x)
        dHdtfun(x, y) = -U * 0.5 * k * cos(k * x)

        errors = map((32, 64, 128, 256)) do n
            h    = lx / n
            grid = StaggeredGrid(Float64, lx, ly, h, h)
            rt   = Runtime(grid)
            topo, mech = TopographicState(grid), MechanicState(grid)

            fill_analytic!(topo.thickness.ice, rt.grid2d, Hfun)
            fill_analytic!(mech.velocity.x_bar, rt.grid2d, (x, y) -> U)
            fill_analytic!(mech.velocity.y_bar, rt.grid2d, (x, y) -> 0.0)

            advect!(topo, mech, UpwindAdvection(), rt)

            exact = analytic_like(topo.thickness.ice_dt, rt.grid2d, dHdtfun)
            maximum(abs, interior(topo.thickness.ice_dt) .- exact)
        end

        rates = convergence_rates(errors)
        @test all(r -> 0.9 < r < 1.1, rates)
    end

    @testset "NoAdvection: ∂H/∂t = ṁ" begin
        grid, rt, topo, mech = setup()

        fill_analytic!(topo.thickness.ice, rt.grid2d, (x, y) -> 500 + 10x)
        fill_analytic!(mech.velocity.x_bar, rt.grid2d, (x, y) -> 3.0)
        fill_analytic!(mech.velocity.y_bar, rt.grid2d, (x, y) -> 3.0)
        fill_analytic!(topo.massbalance.net, rt.grid2d, (x, y) -> 2x - y)

        advect!(topo, mech, NoAdvection(), rt)

        @test all(interior(mech.flux.x) .== 0.0)
        @test all(interior(mech.flux.y) .== 0.0)
        @test interior(topo.thickness.ice_dt) ≈ interior(topo.massbalance.net)
    end

    @testset "scalar mass balance" begin
        grid, rt, topo, mech = setup()
        fill_analytic!(topo.thickness.ice, rt.grid2d, (x, y) -> 500.0)
        fill_analytic!(mech.velocity.x_bar, rt.grid2d, (x, y) -> 0.0)
        fill_analytic!(mech.velocity.y_bar, rt.grid2d, (x, y) -> 0.0)

        mass_flux!(mech.flux.x, mech.flux.y, topo.thickness.ice,
                   mech.velocity.x_bar, mech.velocity.y_bar, CenteredAdvection(), rt)
        thickness_rate!(topo.thickness.ice_dt, mech.flux.x, mech.flux.y, 1.25, rt)

        @test all(interior(topo.thickness.ice_dt) .≈ 1.25)
    end

    # `LevelSetAdvection` is declared but unimplemented: that must be a MethodError, not a
    # silent fallthrough to one of the working schemes.
    @testset "LevelSetAdvection is not implemented" begin
        grid, rt, topo, mech = setup()
        @test_throws MethodError advect!(topo, mech, LevelSetAdvection(), rt)
    end

    # A column grid is the case where `grid2d !== grid`: the fluxes and the thickness live
    # on the shallow grid while the velocity *columns* live on the deep one, so a kernel
    # launched with the wrong launcher would mis-sweep. `x_bar`/`y_bar` are depth-averaged
    # and therefore on `grid2d` in both cases — this pins that the 2D path is selected by
    # the field's grid and not by `nz == 1` accidentally making them the same object.
    @testset "full-column grid: continuity still runs on grid2d" begin
        a, b, U = 100.0, 3.0, 2.0
        grid, rt, topo, mech = setup(; nz = 6)

        @test rt.grid2d !== rt.grid
        @test rt.launch2d !== rt.launch
        @test worksize(rt.launch2d) == (grid.nx + 2, grid.ny + 2, 3)

        fill_analytic!(topo.thickness.ice, rt.grid2d, (x, y) -> a + b * x)
        fill_analytic!(mech.velocity.x_bar, rt.grid2d, (x, y) -> U)
        fill_analytic!(mech.velocity.y_bar, rt.grid2d, (x, y) -> 0.0)

        advect!(topo, mech, CenteredAdvection(), rt)

        @test size(interior(topo.thickness.ice_dt), 3) == 1
        @test all(interior(topo.thickness.ice_dt) .≈ -U * b)
    end

    @testset "depth-integrated grid reuses one launcher" begin
        _, rt, _, _ = setup()
        @test rt.grid2d === rt.grid
        @test rt.launch2d === rt.launch
    end

    @testset "Float32 end to end" begin
        a, b, U = 100.0f0, 3.0f0, 2.0f0
        grid, rt, topo, mech = setup(Float32)

        fill_analytic!(topo.thickness.ice, rt.grid2d, (x, y) -> a + b * x)
        fill_analytic!(mech.velocity.x_bar, rt.grid2d, (x, y) -> U)
        fill_analytic!(mech.velocity.y_bar, rt.grid2d, (x, y) -> 0.0f0)

        advect!(topo, mech, UpwindAdvection(), rt)

        @test eltype(mech.flux.x) === Float32
        @test eltype(topo.thickness.ice_dt) === Float32
        @test eltype(interior(topo.thickness.ice_dt)) === Float32
        @test all(interior(topo.thickness.ice_dt) .≈ -U * b)
    end
end

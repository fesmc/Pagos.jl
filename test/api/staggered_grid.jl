using Pagos
using Test
using KernelAbstractions: CPU

@testset "StaggeredGrid" begin
    @testset "depth-integrated (nz == 1)" begin
        lx, ly, dx, dy = 6000e3, 6000e3, 16e3, 16e3
        grid = StaggeredGrid(Float64, lx, ly, dx, dy)

        @test grid.nx == round(Int, lx / dx)
        @test grid.ny == round(Int, ly / dy)
        @test grid.nz == 1
        @test grid.dx ≈ dx
        @test grid.dy ≈ dy
        @test grid.dz ≈ 1.0                     # placeholder axis, not physically used
        @test size(grid.Lon) == (grid.nx, grid.ny)
        @test size(grid.Lat) == (grid.nx, grid.ny)
        @test size(grid.area) == (grid.nx, grid.ny)
        @test all(grid.area .≈ dx * dy)
        @test all(grid.distortion .≈ 1)
        @test size(grid.basins) == (grid.nx, grid.ny)
        @test size(grid.regions) == (grid.nx, grid.ny)

        # The domain is centred on the origin, but the discretization follows Chmy's
        # C-grid convention (lx/dx cells, centres offset by dx/2 from the domain
        # edges) — deliberately NOT node-for-node compatible with RegularGrid, which
        # places lx/dx + 1 nodes spanning [-lx/2, lx/2] inclusive (see the
        # StaggeredGrid docstring).
        x = collect(grid.x)
        @test x[1] ≈ -lx / 2 + dx / 2
        @test x[end] ≈ lx / 2 - dx / 2
        @test grid.nx == RegularGrid(Float64, lx, ly, dx, dy).nx - 1

        # The size-1 placeholder z-axis is Bounded, NOT Flat: Chmy 0.1.26 exports Flat
        # but defines no methods on it, so it falls through `batch_impl` and makes bc!
        # a MethodError for every field on the grid (see the StaggeredGrid docstring).
        @test connectivity(grid.grid, Dim(3), Side(1)) === Bounded()
        @test connectivity(grid.grid, Dim(3), Side(2)) === Bounded()

        # Pin the reason: bc! must work on a depth-integrated grid.
        f = Field(grid.arch, grid.grid, Center())
        set!(f, 1.0)
        @test (bc!(grid.arch, grid.grid, f => Dirichlet()); true)
    end

    @testset "topology option" begin
        bounded = StaggeredGrid(Float64, 4.0, 4.0, 1.0, 1.0)
        @test connectivity(bounded.grid, Dim(1), Side(1)) === Bounded()
        @test connectivity(bounded.grid, Dim(2), Side(2)) === Bounded()

        connected = StaggeredGrid(Float64, 4.0, 4.0, 1.0, 1.0;
                                  topology = (Connected(), Bounded()))
        @test connectivity(connected.grid, Dim(1), Side(1)) === Connected()
        @test connectivity(connected.grid, Dim(1), Side(2)) === Connected()
        @test connectivity(connected.grid, Dim(2), Side(1)) === Bounded()

        # Periodic and Flat are declared and exported by Chmy 0.1.26 but have no methods
        # anywhere in it, so they fall through dispatch instead of short-circuiting it:
        # `batch_impl` covers Bounded/Connected only, and bc! batches over every axis, so
        # one such axis makes bc! a MethodError for *every* field on the grid. Rejected at
        # construction rather than at the first boundary condition.
        for conn in (Periodic(), Flat())
            @test_throws ErrorException StaggeredGrid(Float64, 4.0, 4.0, 1.0, 1.0;
                                                      topology = (conn, Bounded()))
            @test_throws ErrorException StaggeredGrid(Float64, 4.0, 4.0, 1.0, 1.0;
                                                      topology = (Bounded(), conn))
        end
    end

    @testset "arch constructors agree" begin
        lx, ly, dx, dy = 4.0, 4.0, 1.0, 1.0
        g_default  = StaggeredGrid(Float64, lx, ly, dx, dy)
        g_backend  = StaggeredGrid(CPU(), Float64, lx, ly, dx, dy)
        g_arch     = StaggeredGrid(Arch(CPU()), Float64, lx, ly, dx, dy)

        for g in (g_backend, g_arch)
            @test g.nx == g_default.nx
            @test g.ny == g_default.ny
            @test collect(g.x) ≈ collect(g_default.x)
        end
    end

    @testset "full column (sigma z-axis)" begin
        lx, ly, dx, dy = 6000e3, 6000e3, 16e3, 16e3
        transform = QuadraticSigmaTransform(Float64, 6)
        layering  = CorrectedVerticalLayering(Float64, transform)
        grid      = StaggeredGrid(Float64, lx, ly, dx, dy, layering)

        @test grid.nx == round(Int, lx / dx)
        @test grid.ny == round(Int, ly / dy)
        @test grid.nz == layering.n

        # Chmy derives axis centres from vertices the same way CorrectedVerticalLayering
        # does, so this must match bit-for-bit, not just approximately.
        @test collect(grid.z) == layering.ζ_aa

        @test grid.dx ≈ dx
        @test grid.dy ≈ dy
        @test_throws ErrorException grid.dz     # no single scalar spacing for a sigma axis

        @test connectivity(grid.grid, Dim(3), Side(1)) === Bounded()

        # VerticalLayering does not round-trip through Chmy's vertex-primary axis
        # convention (see the StaggeredGrid docstring) and is deliberately not accepted.
        bad_layering = VerticalLayering(Float64, transform)
        @test_throws MethodError StaggeredGrid(Float64, lx, ly, dx, dy, bad_layering)

        # The sigma axis must be `isbits`: Chmy defines no Adapt rule for StructuredGrid
        # or FunctionAxis, so a grid capturing the layering's Vectors could never be a
        # GPU kernel argument. Pinned on the whole grid, which is what kernels receive.
        @test isbits(axis(grid.grid, Dim(3)))
        @test isbits(grid.grid)
        @test isbits(grid.grid2d)

        # Chmy evaluates the vertex function outside 1:nz+1 — `spacing(ax, Vertex(), k)`
        # is `center(k) - center(k-1)`, so Δz at the bed interface (k == 1) reads
        # vertex(0) and at the surface interface (k == nz+1) reads vertex(nz+2). Those
        # are *interior* points of every z-Vertex field (ε̇_xz, ε̇_yz, w), so indexing
        # ζ_ac directly would be an out-of-range read exactly where DIVA needs a value.
        nz = grid.nz
        for k in 1:(nz + 1)
            @test Δz(grid.grid, Vertex(), 1, 1, k) > 0
        end
        for k in 1:nz
            @test Δz(grid.grid, Center(), 1, 1, k) ≈ layering.ζ_ac[k + 1] - layering.ζ_ac[k]
        end

        # The Launcher sweeps one halo ring (Offset(-1)), so k == 0 and k == nz+2 are
        # reached by ordinary cell-centred kernels too.
        @test zvertex(grid.grid, 0) < layering.ζ_ac[1]
        @test zvertex(grid.grid, nz + 2) > layering.ζ_ac[end]

        # Extrapolating the ghost interfaces must not disturb the interior: cell centres
        # still reproduce ζ_aa bit-for-bit.
        @test collect(grid.z) == layering.ζ_aa
    end

    @testset "full column — Float32 sigma axis" begin
        transform = QuadraticSigmaTransform(Float32, 6)
        layering  = CorrectedVerticalLayering(Float32, transform)
        grid      = StaggeredGrid(Float32, 4f3, 4f3, 1f3, 1f3, layering)

        @test eltype(grid.grid) == Float32
        @test isbits(grid.grid)
        @test collect(grid.z) == layering.ζ_aa
        @test Δz(grid.grid, Vertex(), 1, 1, 1) > 0
        @test Δz(grid.grid, Vertex(), 1, 1, grid.nz + 1) > 0
    end

    @testset "invalid property" begin
        grid = StaggeredGrid(Float64, 4.0, 4.0, 1.0, 1.0)
        @test_throws ErrorException grid.not_a_real_property
    end
end

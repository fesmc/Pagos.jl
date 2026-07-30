using Pagos
using Test

# A field linear in z has an exactly known derivative everywhere, on any axis — no
# discretization error to allow for, so assertions below can use a tight `rtol`/`===`
# rather than a resolution-dependent tolerance.

@testset "∂z_σ" begin
    lx, ly, dx, dy = 4.0, 4.0, 1.0, 1.0
    slope = 3.0

    @testset "exact on a stretched sigma axis, both directions" begin
        layering = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 6))
        grid     = StaggeredGrid(Float64, lx, ly, dx, dy, layering)
        nz       = grid.nz

        # f at Center, linear in zcenter ⟹ ∂z_σ (Center → Vertex) is `slope` at every
        # interface k = 1:nz+1, including the bed (k=1) and surface (k=nz+1) — the
        # indices Chmy's own `∂z` gets wrong except by construction-coincidence.
        f = Field(grid.arch, grid.grid, Center())
        for k in -1:(nz + 2), j in -1:(grid.ny + 2), i in -1:(grid.nx + 2)
            f[i, j, k] = slope * zcenter(grid.grid, k)
        end
        for k in 1:(nz + 1)
            @test ∂z_σ(f, grid.grid, 1, 1, k) ≈ slope rtol = 1e-13
        end

        # f at z-Vertex, linear in zvertex ⟹ ∂z_σ (Vertex → Center) is `slope` at every
        # layer midpoint k = 1:nz. This is the direction where Chmy's own `∂z` is
        # wrong by up to 50% mid-column (see `src/api/sigma_operators.jl`).
        v = Field(grid.arch, grid.grid, (Center(), Center(), Vertex()))
        for k in -1:(nz + 3), j in -1:(grid.ny + 2), i in -1:(grid.nx + 2)
            v[i, j, k] = slope * zvertex(grid.grid, k)
        end
        for k in 1:nz
            @test ∂z_σ(v, grid.grid, 1, 1, k) ≈ slope rtol = 1e-13
        end
    end

    @testset "the bug it fixes is real: Chmy's own ∂z is not exact mid-column" begin
        layering = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 6))
        grid     = StaggeredGrid(Float64, lx, ly, dx, dy, layering)
        f = Field(grid.arch, grid.grid, Center())
        for k in -1:(grid.nz + 2), j in -1:(grid.ny + 2), i in -1:(grid.nx + 2)
            f[i, j, k] = slope * zcenter(grid.grid, k)
        end
        # k = 3 is an interior interface, away from the boundary coincidence at k = 1,
        # nz+1 (see the `∂z_σ` docstring for why those happen to come out right anyway).
        @test !(∂z(f, grid.grid, 1, 1, 3) ≈ slope)
    end

    @testset "reduces to Chmy's own ∂z on a UniformAxis: no regression" begin
        # `UniformAxis` scales by the same precomputed `ax.spacing`/`ax.inv_spacing`
        # regardless of location (`spacing(::UniformAxis, ::Vertex/Center, i)` both
        # return the identical field), so `loc` vs `from` truly cannot matter: exact,
        # bit-for-bit, on the depth-integrated grid's placeholder z-axis.
        grid = StaggeredGrid(Float64, lx, ly, dx, dy)
        f = Field(grid.arch, grid.grid, Center())
        for k in -1:3, j in -1:(grid.ny + 2), i in -1:(grid.nx + 2)
            f[i, j, k] = slope * zcenter(grid.grid, k) + 0.1 * i - 0.2 * j
        end
        @test ∂z_σ(f, grid.grid, 1, 1, 1) === ∂z(f, grid.grid, 1, 1, 1)

        # A *uniformly-spaced* `FunctionAxis` (LinearSigmaTransform) is a subtler case:
        # mathematically the same spacing at Center and Vertex, but `FunctionAxis` has no
        # UniformAxis fast path, so `spacing(ax, Vertex, i) = center(i) - center(i-1)` and
        # `spacing(ax, Center, i) = vertex(i+1) - vertex(i)` take different floating-point
        # paths to the same real number — agreeing to within a few ULP, not bit-for-bit.
        # Still the right regression to pin: it shows `∂z_σ` is a correction, not a
        # behavior change, on the one non-`UniformAxis` case where Chmy already happens to
        # be (numerically) right.
        gridl = StaggeredGrid(Float64, lx, ly, dx, dy,
                               CorrectedVerticalLayering(Float64, LinearSigmaTransform(Float64, 6)))
        nz = gridl.nz
        fl = Field(gridl.arch, gridl.grid, Center())
        for k in -1:(nz + 2), j in -1:(gridl.ny + 2), i in -1:(gridl.nx + 2)
            fl[i, j, k] = slope * zcenter(gridl.grid, k) + 0.1 * i - 0.2 * j
        end
        for k in 1:(nz + 1)
            @test ∂z_σ(fl, gridl.grid, 1, 1, k) ≈ ∂z(fl, gridl.grid, 1, 1, k) rtol = 1e-13
        end
    end

    @testset "zero-allocation" begin
        layering = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 6))
        grid     = StaggeredGrid(Float64, lx, ly, dx, dy, layering)
        f = Field(grid.arch, grid.grid, Center())
        set!(f, 1.0)
        ∂z_σ(f, grid.grid, 1, 1, 2)   # compile
        @test (@allocated ∂z_σ(f, grid.grid, 1, 1, 2)) == 0
    end
end

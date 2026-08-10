using Pagos
using Test
using Random

# `fill_from_grid!`/`fill_from_grid3d!`: a halo-aware bulk fill for seeding a field
# directly from an external (small, host-side) array — the x/y halo ring is filled by
# clamping to the nearest interior column/row, matching what a naive
# `f[i, j, k] = data[clamp(i, 1, ni), clamp(j, 1, nj)]` loop over the halo would produce,
# but built from bulk array ops (`copyto!`) so it works on a device field without scalar
# indexing. See `pagos-roadmap/memreduce.md`, "Build device states without a host
# duplicate".
#
# Correctness is checked against that scalar loop directly, generalized to whatever halo
# the field actually has (`Field`'s `halo` type parameter H gives `2H` reachable ghost
# rings per side — see `Chmy.Fields.Field`), rather than the docs example's original
# hardcoded halo=1 loop.

_ref_halo_xy(f) = (h = halo(f); h isa Tuple ? (h[1], h[2]) : (h, h))

function ref_fill_from_grid!(f, data)
    ni, nj = size(data)
    hx, hy = _ref_halo_xy(f)
    for k in axes(interior(f), 3), j in (1-2hy):(nj+2hy), i in (1-2hx):(ni+2hx)
        ic, jc = clamp(i, 1, ni), clamp(j, 1, nj)
        f[i, j, k] = data[ic, jc]
    end
    return f
end

function ref_fill_from_grid3d!(f, data)
    ni, nj, nk = size(data)
    hx, hy = _ref_halo_xy(f)
    for k = 1:nk, j in (1-2hy):(nj+2hy), i in (1-2hx):(ni+2hx)
        ic, jc = clamp(i, 1, ni), clamp(j, 1, nj)
        f[i, j, k] = data[ic, jc, k]
    end
    return f
end

@testset "fill_from_grid!/fill_from_grid3d!" begin
    Random.seed!(42)

    @testset "matches scalar-loop reference (scalar halo=1, grid2d)" begin
        grid = StaggeredGrid(Float64, 40.0, 30.0, 1.0, 1.0)
        ni, nj = grid.nx, grid.ny
        data = rand(ni, nj)

        f_new = Field(grid.arch, grid.grid2d, (Center(), Center(), Center()), Float64; halo = 1)
        f_ref = Field(grid.arch, grid.grid2d, (Center(), Center(), Center()), Float64; halo = 1)

        fill_from_grid!(f_new, data)
        ref_fill_from_grid!(f_ref, data)

        @test parent(f_new) == parent(f_ref)
    end

    @testset "matches reference on a non-square grid, halo=2" begin
        grid = StaggeredGrid(Float64, 25.0, 60.0, 1.0, 1.0)
        ni, nj = grid.nx, grid.ny
        data = rand(ni, nj)

        f_new = Field(grid.arch, grid.grid2d, (Center(), Center(), Center()), Float64; halo = 2)
        f_ref = Field(grid.arch, grid.grid2d, (Center(), Center(), Center()), Float64; halo = 2)

        fill_from_grid!(f_new, data)
        ref_fill_from_grid!(f_ref, data)

        @test parent(f_new) == parent(f_ref)
    end

    @testset "matches reference with per-axis halo (hx, hy, 0), item-4 style" begin
        grid = StaggeredGrid(Float64, 25.0, 60.0, 1.0, 1.0)
        ni, nj = grid.nx, grid.ny
        data = rand(ni, nj)

        f_new =
            Field(grid.arch, grid.grid2d, (Center(), Center(), Center()), Float64; halo = (1, 1, 0))
        f_ref =
            Field(grid.arch, grid.grid2d, (Center(), Center(), Center()), Float64; halo = (1, 1, 0))

        fill_from_grid!(f_new, data)
        ref_fill_from_grid!(f_ref, data)

        @test parent(f_new) == parent(f_ref)
    end

    @testset "matches reference with unequal per-axis halo (hx != hy)" begin
        grid = StaggeredGrid(Float64, 25.0, 60.0, 1.0, 1.0)
        ni, nj = grid.nx, grid.ny
        data = rand(ni, nj)

        f_new =
            Field(grid.arch, grid.grid2d, (Center(), Center(), Center()), Float64; halo = (2, 1, 0))
        f_ref =
            Field(grid.arch, grid.grid2d, (Center(), Center(), Center()), Float64; halo = (2, 1, 0))

        fill_from_grid!(f_new, data)
        ref_fill_from_grid!(f_ref, data)

        @test parent(f_new) == parent(f_ref)
    end

    @testset "fill_from_grid3d! matches scalar-loop reference" begin
        nk = 6
        layering = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, nk))
        grid = StaggeredGrid(Float64, 25.0, 30.0, 1.0, 1.0, layering)
        ni, nj = grid.nx, grid.ny
        data = rand(ni, nj, nk)

        f_new = Field(grid.arch, grid.grid, (Center(), Center(), Center()), Float64; halo = 1)
        f_ref = Field(grid.arch, grid.grid, (Center(), Center(), Center()), Float64; halo = 1)

        fill_from_grid3d!(f_new, data)
        ref_fill_from_grid3d!(f_ref, data)

        @test parent(f_new) == parent(f_ref)
    end

    @testset "leaves the z halo untouched (scalar halo case)" begin
        grid = StaggeredGrid(Float64, 20.0, 20.0, 1.0, 1.0)
        ni, nj = grid.nx, grid.ny
        data = rand(ni, nj)
        f = Field(grid.arch, grid.grid2d, (Center(), Center(), Center()), Float64; halo = 1)
        fill_from_grid!(f, data)
        @test all(==(0.0), parent(f)[:, :, 1])
        @test all(==(0.0), parent(f)[:, :, 2])
    end

    @testset "fill_from_grid! on a real MechanicState field (grid2d, item-4 shrunk halo)" begin
        grid = StaggeredGrid(Float64, 25.0, 30.0, 1.0, 1.0)
        mech = MechanicState(grid)
        ni, nj = grid.nx, grid.ny
        data = rand(ni, nj)

        fill_from_grid!(mech.friction.beta, data)
        @test dropdims(asarray(mech.friction.beta); dims = 3) == data
    end
end

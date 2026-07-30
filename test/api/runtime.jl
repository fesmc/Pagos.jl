using Pagos
using Test
using KernelAbstractions: CPU, @kernel, @index

@kernel inbounds = true function _fill_x!(f, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    f[I...] = xcenters(grid)[I[1]]
end

@kernel inbounds = true function _fill_const!(f, v, O)
    I = @index(Global, NTuple)
    I = I + O
    f[I...] = v
end

# A Chmy-idiomatic kernel: the Chmy grid is an ordinary kernel argument, and a Chmy
# operator dispatches on it. Pins that `rt.grid` is directly usable this way.
@kernel inbounds = true function _ddx!(dfdx, f, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    dfdx[I...] = ∂x(f, grid, I...)
end

@kernel inbounds = true function _fill_z!(f, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    f[I...] = zcoord(grid, location(f), I...)
end

@kernel inbounds = true function _ddz!(dfdz, f, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    dfdz[I...] = ∂z(f, grid, I...)
end

@testset "Runtime" begin
    lx, ly, dx, dy = 8.0, 8.0, 1.0, 1.0

    @testset "construction" begin
        sgrid = StaggeredGrid(Float64, lx, ly, dx, dy)
        rt    = Runtime(sgrid)

        # `grid` is the bare Chmy grid, NOT the StaggeredGrid wrapper: Chmy's
        # operators, Launcher and bc! all dispatch on Chmy.StructuredGrid, and the
        # wrapper has no Adapt rule so it could never enter a kernel.
        @test rt.grid === sgrid.grid
        @test rt.grid isa Chmy.Grids.StructuredGrid
        @test rt.arch === sgrid.arch
        @test rt.launch isa Launcher
        @test Pagos.KernelAbstractions.get_backend(rt) === CPU()

        # Chmy's Launcher sweeps interior + one halo ring, in all three dims (the
        # z-axis is the size-1 placeholder here). Pinned because Phase 3/4 kernels
        # depend on this worksize convention.
        @test worksize(rt.launch) == (sgrid.nx + 2, sgrid.ny + 2, sgrid.nz + 2)
        @test outer_width(rt.launch) === nothing
    end

    @testset "outer_width forwarded" begin
        rt = Runtime(StaggeredGrid(Float64, lx, ly, dx, dy); outer_width = (2, 2, 1))
        @test outer_width(rt.launch) == (2, 2, 1)
    end

    @testset "launches with Chmy's documented signature" begin
        sgrid = StaggeredGrid(Float64, lx, ly, dx, dy)
        rt    = Runtime(sgrid)
        f     = Field(rt.arch, rt.grid, Center())

        rt.launch(rt.arch, rt.grid, _fill_x! => (f, rt.grid))

        # interior holds the cell-centre x-coordinate of each point
        @test interior(f) ≈ repeat(collect(sgrid.x), 1, sgrid.ny, sgrid.nz)

        # ...and the one halo ring the Launcher sweeps was written too, rather than
        # left at zero. Note this is `interior(f; with_halo=true)`, not `parent(f)`:
        # Chmy allocates `size .+ 4 .* halo`, i.e. TWO ghost rings per side for the
        # default halo=1, while the Launcher's worksize covers only the inner one.
        @test all(interior(f; with_halo = true) .!= 0)
        @test count(iszero, parent(f)) > 0      # outer ring: allocated, never swept
    end

    @testset "rt.grid drives Chmy operators inside a kernel" begin
        sgrid = StaggeredGrid(Float64, lx, ly, dx, dy)
        rt    = Runtime(sgrid)
        f     = Field(rt.arch, rt.grid, Center())
        dfdx  = Field(rt.arch, rt.grid, (Vertex(), Center(), Center()))

        set!(f, rt.grid, (x, y, z) -> 2x)
        rt.launch(rt.arch, rt.grid, _ddx! => (dfdx, f, rt.grid))

        # ∂x of a field varying as 2x is 2 everywhere the stencil is well-posed
        # (interior vertices; the outermost ones read an unset ghost).
        @test all(interior(dfdx)[2:(end - 1), :, :] .≈ 2.0)
    end

    @testset "launch applies boundary conditions" begin
        rt = Runtime(StaggeredGrid(Float64, lx, ly, dx, dy))
        f  = Field(rt.arch, rt.grid, Center())

        # Fill with a constant, then let a Dirichlet(0) bc fill the ghost cells.
        rt.launch(rt.arch, rt.grid, _fill_const! => (f, 1.0);
                  bc = batch(rt.grid, f => Dirichlet()))

        @test all(interior(f) .≈ 1.0)

        # Dirichlet(0) on a Center field is enforced halfway between the last interior
        # point and its ghost, so the two must average to zero. Corners are excluded:
        # bc! fills each face's ghosts from interior points only, so the cells where
        # two ghost rings meet keep whatever the kernel's halo sweep left there.
        halo_x = interior(f; with_halo = true)[1, 2:(end - 1), 2:(end - 1)]
        @test all(halo_x .≈ -1.0)
    end

    # Every other Runtime test runs on a depth-integrated grid, whose z-axis is a plain
    # UniformAxis. The column grid's sigma FunctionAxis is the one Chmy evaluates outside
    # 1:nz+1 — both because the Launcher sweeps a halo ring (Offset(-1)) and because
    # `spacing(ax, Vertex(), k)` reaches vertex(k-1). A kernel reading z geometry is what
    # catches that; a host-side `collect(grid.z)` never does.
    @testset "launches on a full-column grid" begin
        layering = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 6))
        sgrid    = StaggeredGrid(Float64, lx, ly, dx, dy, layering)
        rt       = Runtime(sgrid)

        f = Field(rt.arch, rt.grid, Center())
        rt.launch(rt.arch, rt.grid, _fill_z! => (f, rt.grid))
        @test interior(f)[1, 1, :] ≈ layering.ζ_aa

        # ∂z of a z-Center field lands on z-Vertex, whose first and last layers are the
        # bed and surface interfaces — exactly where ε̇_xz/ε̇_yz live, and exactly the
        # indices that reach vertex(0) / vertex(nz+2). Every interface must be finite.
        dfdz = Field(rt.arch, rt.grid, (Center(), Center(), Vertex()))
        rt.launch(rt.arch, rt.grid, _ddz! => (dfdz, f, rt.grid))
        @test all(isfinite, interior(dfdz))
        @test size(interior(dfdz), 3) == sgrid.nz + 1

        # NOTE: this is *not* ≈ 1 even though f == ζ. Chmy's `∂` differences across the
        # destination stencil but scales by `iΔ(grid, location(f), ...)` — the spacing at
        # the *source* location. For a Center→Vertex ∂z at k that is ζ_ac[k+1] - ζ_ac[k],
        # not the ζ_aa[k] - ζ_aa[k-1] the difference actually spans. The two coincide on a
        # uniform axis, which is why nothing else here notices; on the stretched sigma
        # axis they do not. This is the bug `∂z_σ` (`src/api/sigma_operators.jl`) works
        # around; pinned here, on Chmy's raw `∂z`, so a Chmy bump that changes the
        # convention is visible (at which point `∂z_σ` may become redundant).
        k = 3
        @test interior(dfdz)[1, 1, k] ≈
              (layering.ζ_aa[k] - layering.ζ_aa[k - 1]) / (layering.ζ_ac[k + 1] - layering.ζ_ac[k])
    end

    @testset "Float32 grid" begin
        sgrid = StaggeredGrid(Float32, lx, ly, dx, dy)
        rt    = Runtime(sgrid)
        f     = Field(rt.arch, rt.grid, Center())

        @test eltype(f) === Float32
        rt.launch(rt.arch, rt.grid, _fill_x! => (f, rt.grid))
        @test interior(f) ≈ repeat(collect(sgrid.x), 1, sgrid.ny, sgrid.nz)
        @test eltype(interior(f)) === Float32
    end
end

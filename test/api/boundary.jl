using Pagos
using Test

# The plain-array API boundary: `asarray` out, `setdata!` in, both working on either
# representation (Chmy `Field`s or plain arrays) so that public/user-facing code never has
# to name a Chmy type. See `pagos-roadmap/chmy.md`, Phase 2.

@testset "plain-array API boundary" begin
    lx, ly, dx, dy = 8.0, 8.0, 1.0, 1.0
    sgrid = StaggeredGrid(Float64, lx, ly, dx, dy)
    rgrid = RegularGrid(Float64, lx, ly, dx, dy)

    @testset "asarray on both representations" begin
        f = Field(sgrid.arch, sgrid.grid, Center())
        A = zeros(3, 4)

        @test asarray(f) === interior(f)
        @test size(asarray(f)) == (sgrid.nx, sgrid.ny, sgrid.nz)
        @test asarray(A) === A          # a no-op on plain arrays, not a copy

        # It is a view: writing through it lands in the field, and the halo is untouched.
        fill!(asarray(f), 7.0)
        @test all(interior(f) .≈ 7.0)
        @test all(iszero, parent(f)[1, :, :])
    end

    # A Chmy Field IS an AbstractArray (`AbstractField{T,N,L} <: AbstractArray{T,N}`), so
    # the array method would also apply to it — the field method must be the one selected,
    # or `asarray` would hand back the field itself and the boundary would leak.
    @testset "asarray dispatches to the field method, not the array fallback" begin
        f = Field(sgrid.arch, sgrid.grid, Center())
        @test f isa AbstractArray
        @test !(asarray(f) isa Chmy.Fields.AbstractField)
        @test asarray(f) isa SubArray
    end

    @testset "asarray over a whole state, either grid kind" begin
        for grid in (sgrid, rgrid)
            topo = TopographicState(grid)
            mech = MechanicState(grid)
            for f in (topo.thickness.ice, topo.elevation.surface_dx,
                      mech.flux.x, mech.velocity.x, mech.strainrate.xz)
                @test asarray(f) isa AbstractArray
                @test !(asarray(f) isa Chmy.Fields.AbstractField)
            end
        end
    end

    @testset "setdata! from a scalar" begin
        for grid in (sgrid, rgrid)
            topo = TopographicState(grid)
            @test setdata!(topo.thickness.ice, 1500.0) === topo.thickness.ice
            @test all(asarray(topo.thickness.ice) .≈ 1500.0)
        end
    end

    @testset "setdata! from a matching array" begin
        topo = TopographicState(sgrid)
        A = reshape(collect(1.0:(sgrid.nx * sgrid.ny)), sgrid.nx, sgrid.ny, 1)
        setdata!(topo.thickness.ice, A)
        @test asarray(topo.thickness.ice) == A
    end

    # The common ergonomic case: the user has an (nx, ny) matrix and the depth-integrated
    # field is (nx, ny, 1). Trailing singletons carry no data, so this is accepted.
    @testset "setdata! accepts a 2D matrix for a depth-integrated field" begin
        topo = TopographicState(sgrid)
        M = reshape(collect(1.0:(sgrid.nx * sgrid.ny)), sgrid.nx, sgrid.ny)
        setdata!(topo.thickness.ice, M)
        @test asarray(topo.thickness.ice)[:, :, 1] == M

        # ...and the reverse pairing, writing a 3-D source into a 2-D destination.
        dst = zeros(sgrid.nx, sgrid.ny)
        setdata!(dst, reshape(M, sgrid.nx, sgrid.ny, 1))
        @test dst == M
    end

    # The reason `setdata!` exists rather than calling `Chmy.set!` directly: Chmy's
    # `set!(f, A::AbstractArray)` is a bare `copyto!`, which copies by *linear* index and
    # only errors when the source is too large. Each case below is silent data corruption
    # through `Chmy.set!` and a `DimensionMismatch` through `setdata!`.
    @testset "setdata! rejects what Chmy.set! accepts silently" begin
        f = Field(sgrid.arch, sgrid.grid, Center())

        @testset "too small: Chmy.set! partial-fills without error" begin
            set!(f, 0.0)
            Chmy.set!(f, ones(3, 3))                        # no error raised
            @test sum(interior(f)) == 9.0                   # only 9 of 64 cells written
            @test count(iszero, interior(f)) == 55          # ...the rest silently stale
            @test_throws DimensionMismatch setdata!(f, ones(3, 3))
        end

        @testset "right length, wrong shape: Chmy.set! reinterprets linearly" begin
            B = reshape(collect(1.0:(sgrid.nx * sgrid.ny)), 4, 16)
            Chmy.set!(f, B)                                 # no error raised
            @test interior(f)[1, 2, 1] == 9.0               # read from B[1, 3], linearly
            @test_throws DimensionMismatch setdata!(f, B)
        end

        # A `RegularGrid` of the same (lx, ly, dx, dy) is one node larger per axis, so data
        # carried over from the old grid is off by one in each direction. Chmy does at
        # least error here (the source is too large), but with a BoundsError that names no
        # sizes; the checked path says what does not fit.
        @testset "RegularGrid-sized data: off by one per horizontal axis" begin
            @test rgrid.nx == sgrid.nx + 1
            legacy = ones(rgrid.nx, rgrid.ny)
            @test_throws BoundsError Chmy.set!(f, legacy)
            @test_throws DimensionMismatch setdata!(f, legacy)
        end
    end

    @testset "setdata! converts element type and preserves the field's" begin
        grid32 = StaggeredGrid(Float32, lx, ly, dx, dy)
        topo   = TopographicState(grid32)

        setdata!(topo.thickness.ice, 1234.5)                       # a Float64 scalar
        @test eltype(asarray(topo.thickness.ice)) === Float32
        @test asarray(topo.thickness.ice)[1, 1, 1] ≈ 1234.5f0

        setdata!(topo.thickness.ice, ones(grid32.nx, grid32.ny))   # a Float64 array
        @test eltype(asarray(topo.thickness.ice)) === Float32
        @test all(asarray(topo.thickness.ice) .≈ 1.0f0)
    end

    @testset "setdata! from another field copies interior to interior" begin
        src = Field(sgrid.arch, sgrid.grid, Center())
        dst = Field(sgrid.arch, sgrid.grid, Center())
        set!(src, 3.0)
        parent(src) .= ifelse.(parent(src) .== 3.0, 3.0, -99.0)     # poison the halo

        setdata!(dst, src)
        @test all(interior(dst) .≈ 3.0)
        @test !any(interior(dst) .≈ -99.0)      # the halo was not dragged along
    end

    # A column field is where the third dimension is real, so the singleton-squeeze must
    # not paper over a layer-count mismatch.
    @testset "column fields: the vertical dimension is checked" begin
        layering = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 6))
        grid = StaggeredGrid(Float64, lx, ly, dx, dy, layering)
        mech = MechanicState(grid)
        η    = mech.material.viscosity                # aa, column

        @test size(asarray(η)) == (grid.nx, grid.ny, grid.nz)
        setdata!(η, ones(grid.nx, grid.ny, grid.nz))
        @test all(asarray(η) .≈ 1.0)
        @test_throws DimensionMismatch setdata!(η, ones(grid.nx, grid.ny, grid.nz - 1))
        @test_throws DimensionMismatch setdata!(η, ones(grid.nx, grid.ny))

        # A z-Vertex field has nz + 1 layers, not nz — the asymmetry the state docs warn
        # about. It must not silently accept nz layers of data.
        @test size(asarray(mech.strainrate.xz), 3) == grid.nz + 1
        @test_throws DimensionMismatch setdata!(mech.strainrate.xz,
                                               ones(grid.nx + 1, grid.ny, grid.nz))
    end

    # The one existing diagnostic that took an AbstractArray and is now field-aware: a
    # Field IS an AbstractArray, so it already "worked", but by scalar-indexing the field
    # instead of reducing over the underlying array, and it must look at the interior only.
    @testset "hasnan goes through asarray" begin
        f = Field(sgrid.arch, sgrid.grid, Center())
        set!(f, 0.0)
        @test !Pagos.hasnan(f)

        parent(f) .= ifelse.(parent(f) .== 0.0, 0.0, NaN)
        parent(f)[1, 1, 1] = NaN                # an outer-halo cell, never swept
        @test !Pagos.hasnan(f)                  # halo NaNs are not solution NaNs

        interior(f)[2, 2, 1] = NaN
        @test Pagos.hasnan(f)
        @test Pagos.hasnan([1.0, NaN])          # plain arrays unchanged
        @test !Pagos.hasnan([1.0, 2.0])
    end

    # The round trip a user of the public API actually performs: plain arrays in, plain
    # arrays out, with no Chmy type named at any point.
    @testset "array -> state -> array round trip" begin
        topo = TopographicState(sgrid)
        H_in = 500 .+ 100 .* rand(sgrid.nx, sgrid.ny)

        setdata!(topo.thickness.ice, H_in)
        H_out = asarray(topo.thickness.ice)

        @test H_out[:, :, 1] == H_in
        @test H_out isa AbstractArray{Float64, 3}

        # ...and the identical call sequence works on the plain-array representation,
        # sized against its own grid (`RegularGrid` is one node larger per axis).
        topo_plain = TopographicState(rgrid)
        H_legacy   = 500 .+ 100 .* rand(rgrid.nx, rgrid.ny)
        setdata!(topo_plain.thickness.ice, H_legacy)
        @test asarray(topo_plain.thickness.ice) == H_legacy
    end
end

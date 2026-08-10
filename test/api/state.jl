using Pagos
using Test
using KernelAbstractions: @kernel, @index

# The whole point of location-typed states is that a wrong location is a *type* error
# rather than a silent half-cell shift, so these tests pin the location of every field
# group, and then check that the layout is consistent with the operators that will read
# it: `∂x` of a field at `acx` must land exactly on the `aa` field it is written into.
@kernel inbounds = true function _ddx!(dfdx, f, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    dfdx[I...] = ∂x(f, grid, I...)
end

@kernel inbounds = true function _ddy!(dfdy, f, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    dfdy[I...] = ∂y(f, grid, I...)
end

@testset "state" begin
    lx, ly, dx, dy = 8.0, 6.0, 1.0, 1.0
    layering = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, 6))

    @testset "field locations" begin
        grid = StaggeredGrid(Float64, lx, ly, dx, dy, layering)
        topo = TopographicState(grid)
        mech = MechanicState(grid)
        therm = ThermodynamicState(grid)
        mat = MaterialState(grid)

        # Topography: everything at aa, except the surface gradients on the faces.
        @test location(topo.thickness.ice) == Pagos.NODE_AA
        @test location(topo.mask.is_ice) == Pagos.NODE_AA
        @test location(topo.elevation.surface) == Pagos.NODE_AA
        @test location(topo.elevation.surface_dx) == Pagos.NODE_ACX
        @test location(topo.elevation.surface_dy) == Pagos.NODE_ACY
        @test eltype(parent(topo.mask.is_grounded)) === Bool

        # Velocity: u on acx, v on acy, w on the layer interfaces at aa.
        @test location(mech.velocity.x) == Pagos.NODE_ACX
        @test location(mech.velocity.y) == Pagos.NODE_ACY
        @test location(mech.velocity.z) == Pagos.NODE_AA_AC
        # ... and each velocity gradient sits where the operator puts it.
        @test location(mech.velocity.x_dx) == Pagos.NODE_AA
        @test location(mech.velocity.x_dy) == Pagos.NODE_AB
        @test location(mech.velocity.x_dz) == Pagos.NODE_ACX_AC
        @test location(mech.velocity.y_dx) == Pagos.NODE_AB
        @test location(mech.velocity.y_dy) == Pagos.NODE_AA
        @test location(mech.velocity.y_dz) == Pagos.NODE_ACY_AC
        @test location(mech.velocity.z_dx) == Pagos.NODE_ACX_AC
        @test location(mech.velocity.z_dy) == Pagos.NODE_ACY_AC
        @test location(mech.velocity.z_dz) == Pagos.NODE_AA

        # Tensors: normal components at aa, in-plane shear on the corners, vertical
        # shear on the layer interfaces at the corresponding face. Same convention as
        # Chmy's own `TensorField{3}`.
        for tensor in (mech.strainrate, mech.stress)
            @test location(tensor.xx) == Pagos.NODE_AA
            @test location(tensor.yy) == Pagos.NODE_AA
            @test location(tensor.zz) == Pagos.NODE_AA
            @test location(tensor.xy) == Pagos.NODE_AB
            @test location(tensor.yx) == Pagos.NODE_AB
            @test location(tensor.xz) == Pagos.NODE_ACX_AC
            @test location(tensor.zx) == Pagos.NODE_ACX_AC
            @test location(tensor.yz) == Pagos.NODE_ACY_AC
            @test location(tensor.zy) == Pagos.NODE_ACY_AC
            @test location(tensor.effective) == Pagos.NODE_AA
        end

        # Driving and basal stress live on the velocity faces they accelerate.
        @test location(mech.stress.driving_x) == Pagos.NODE_ACX
        @test location(mech.stress.driving_y) == Pagos.NODE_ACY
        @test location(mech.stress.base_x) == Pagos.NODE_ACX
        @test location(mech.stress.base_y) == Pagos.NODE_ACY

        # Depth-averaged velocity and its gradients follow the same algebra as the
        # column ones, on the single-layer grid.
        @test location(mech.velocity.depthaverage_x) == Pagos.NODE_ACX
        @test location(mech.velocity.depthaverage_x_dy) == Pagos.NODE_AB
        @test location(mech.velocity.depthaverage_y_dx) == Pagos.NODE_AB
        @test location(mech.friction.beta) == Pagos.NODE_AA

        @test location(mat.eta_ice) == Pagos.NODE_AA
        @test location(therm.temperature.ice) == Pagos.NODE_AA
    end

    @testset "column vs depth-integrated dimensioning" begin
        grid = StaggeredGrid(Float64, lx, ly, dx, dy, layering)
        (; nx, ny, nz) = grid
        @test (nx, ny, nz) == (8, 6, 6)

        mech = MechanicState(grid)
        mat = MaterialState(grid)
        therm = ThermodynamicState(grid)
        topo = TopographicState(grid)

        # Column fields carry all nz layers; a Vertex z-location carries nz + 1.
        @test size(mech.velocity.x) == (nx + 1, ny, nz)
        @test size(mech.velocity.y) == (nx, ny + 1, nz)
        @test size(mech.velocity.z) == (nx, ny, nz + 1)
        @test size(mech.strainrate.xy) == (nx + 1, ny + 1, nz)
        @test size(mech.strainrate.xz) == (nx + 1, ny, nz + 1)
        @test size(mat.eta_ice) == (nx, ny, nz)
        @test size(therm.temperature.ice) == (nx, ny, nz)

        # Depth-integrated fields are single-layer, on `grid2d` — but still 3D, so they
        # are indexed `f[i, j, 1]`.
        @test size(mech.stress.driving_x) == (nx + 1, ny, 1)
        @test size(mech.velocity.depthaverage_x) == (nx + 1, ny, 1)
        @test size(mech.friction.beta) == (nx, ny, 1)
        @test size(topo.thickness.ice) == (nx, ny, 1)
        @test size(therm.temperature.ice_surface) == (nx, ny, 1)
        @test size(mat.eta_depth_averaged) == (nx, ny, 1)
        @test ndims(topo.thickness.ice) == 3
    end

    @testset "depth-integrated grid" begin
        grid = StaggeredGrid(Float64, lx, ly, dx, dy)
        (; nx, ny, nz) = grid
        @test nz == 1
        # With a single layer there is nothing to separate: both grids are the same one.
        @test grid.grid2d === grid.grid

        mech = MechanicState(grid)
        @test size(mech.velocity.x) == (nx + 1, ny, 1)
        @test size(mech.stress.driving_x) == (nx + 1, ny, 1)
        # The layer-interface fields still have both interfaces (bed and surface) even
        # for a single layer — vertical shear is not a cell-centred quantity.
        @test size(mech.velocity.z) == (nx, ny, 2)
        @test size(mech.strainrate.xz) == (nx + 1, ny, 2)

        column = StaggeredGrid(Float64, lx, ly, dx, dy, layering)
        @test column.grid2d !== column.grid
        @test size(column.grid2d, Center())[1:2] == size(column.grid, Center())[1:2]
        @test size(column.grid2d, Center())[3] == 1
        @test xcenters(column.grid2d) == xcenters(column.grid)
    end

    @testset "layout is operator-consistent" begin
        # Not a restatement of the locations: this checks that Chmy's operators map one
        # state field onto another. ∂x of `velocity.x` (acx) must land exactly on
        # `velocity.x_dx` (aa), and ∂y of it on `velocity.x_dy` (ab) — if either
        # location were wrong the field written into would be the wrong size, or the
        # values would be shifted by half a cell.
        grid = StaggeredGrid(Float64, lx, ly, dx, dy, layering)
        (; ny) = grid
        rt = Runtime(grid)
        mech = MechanicState(grid)
        (; velocity) = mech

        set!(velocity.x, grid.grid, (x, y, z) -> 2x + 3y)
        rt.launch(rt.arch, rt.grid, _ddx! => (velocity.x_dx, velocity.x, rt.grid))
        rt.launch(rt.arch, rt.grid, _ddy! => (velocity.x_dy, velocity.x, rt.grid))

        # ∂x maps acx → aa: both cells feeding a centre are interior, so the whole
        # field is exact without any halo fill.
        @test all(≈(2), interior(velocity.x_dx))
        # ∂y maps acx → ab, and the outermost y-vertices read a y-halo cell of `u` that
        # nothing has filled — so the layout is exact on the vertices that have interior
        # neighbours on both sides, and the boundary rows are a Phase 4 (`bc!`) question,
        # not a layout one.
        @test all(≈(3), view(interior(velocity.x_dy), :, 2:ny, :))
        @test !all(≈(3), interior(velocity.x_dy))
    end

    @testset "eltype and halo follow the grid" begin
        grid32 = StaggeredGrid(Float32, lx, ly, dx, dy, layering)
        mech32 = MechanicState(grid32)
        topo32 = TopographicState(grid32)
        @test eltype(parent(mech32.velocity.x)) === Float32
        @test eltype(parent(mech32.stress.driving_x)) === Float32
        @test eltype(parent(topo32.elevation.surface)) === Float32
        @test eltype(parent(topo32.mask.is_ice)) === Bool     # masks stay boolean

        grid = StaggeredGrid(Float64, lx, ly, dx, dy, layering)
        # Column fields (on `grid.grid`): halo is the scalar the caller passed, on every axis.
        @test Chmy.Fields.halo(MechanicState(grid).velocity.x) == 1
        @test Chmy.Fields.halo(MechanicState(grid; halo = 2).velocity.x) == 2
        @test Chmy.Fields.halo(ThermodynamicState(grid; halo = 2).enthalpy.ice) == 2
        @test Chmy.Fields.halo(MaterialState(grid; halo = 2).eta_ice) == 2

        # Depth-integrated fields (on `grid.grid2d`): `z` gets none — `_halo2d`, and
        # `pagos-roadmap/memreduce.md`'s per-axis-halo item. `TopographicState` is entirely
        # depth-integrated, so `thickness.ice` here stands for the whole struct.
        @test Chmy.Fields.halo(MechanicState(grid).friction.beta) == (1, 1, 0)
        @test Chmy.Fields.halo(MechanicState(grid; halo = 2).friction.beta) == (2, 2, 0)
        @test Chmy.Fields.halo(TopographicState(grid; halo = 2).thickness.ice) == (2, 2, 0)
        @test Chmy.Fields.halo(ThermodynamicState(grid; halo = 2).heatflux_base_ice) ==
              (2, 2, 0)
        @test Chmy.Fields.halo(MaterialState(grid; halo = 2).eta_depth_averaged) == (2, 2, 0)

        # The exception: `MechanicState`'s `ACX2`/`ACY2`-typed fields keep a full,
        # every-axis halo, because Chmy's `bc!` sweeps every axis *other* than the one it
        # fills at the grid's own `size .+ 2` — so a `(h, h, 0)`-halo field goes out of
        # bounds even when only `x`/`y` are being filled, and `velocity.depthaverage_x`/`y`
        # are `bc!`'d every PT iteration (`pseudo_transient!`). See the note on
        # `MechanicState`'s constructor. `flux.x`/`stress.driving_x`/`velocity.base_x` share
        # the *same* type parameter as `depthaverage_x` and so share its halo too, whether
        # or not each one is itself ever `bc!`'d.
        @test Chmy.Fields.halo(MechanicState(grid).velocity.depthaverage_x) == 1
        @test Chmy.Fields.halo(MechanicState(grid; halo = 2).velocity.depthaverage_x) == 2
        @test Chmy.Fields.halo(MechanicState(grid).flux.x) == 1
        @test Chmy.Fields.halo(MechanicState(grid).stress.driving_x) == 1
        @test Chmy.Fields.halo(MechanicState(grid).velocity.base_x) == 1
        @test Chmy.Fields.halo(MechanicState(grid).velocity.surface_x) == 1
        @test Chmy.Fields.halo(MechanicState(grid).velocity.depthaverage_x_dz) == 1
    end

    @testset "adapt" begin
        # `Field` carries its own Adapt rule, so a state of Fields stays adaptable —
        # this is what a GPU launch relies on.
        grid = StaggeredGrid(Float64, lx, ly, dx, dy, layering)
        mech = Pagos.Adapt.adapt(Array, MechanicState(grid))
        @test mech isa MechanicState
        @test location(mech.velocity.x) == Pagos.NODE_ACX
        @test size(mech.velocity.x) == size(MechanicState(grid).velocity.x)
    end

    @testset "RegularGrid states are unchanged" begin
        # The `Field` constructors are additive: the plain-array path keeps its own
        # shapes (2D matrices for the depth-integrated fields, no halo, no location).
        grid = RegularGrid(Float64, lx, ly, dx, dy)
        (; nx, ny, nz) = grid
        mech = MechanicState(grid)
        topo = TopographicState(grid)
        therm = ThermodynamicState(grid)
        mat = MaterialState(grid)

        @test mech.velocity.x isa Array{Float64, 3}
        @test size(mech.velocity.x) == (nx, ny, nz)
        @test size(mech.stress.driving_x) == (nx, ny)
        @test size(topo.thickness.ice) == (nx, ny)
        @test size(topo.elevation.surface_dx) == (nx, ny)
        @test size(therm.temperature.ice) == (nx, ny)
        @test size(mat.eta_ice) == (nx, ny, nz)
    end
end

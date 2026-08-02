using Pagos
using Test

include("../test_helpers/chmy.jl")

# Ice masks: derivation, the node-class reduction rules, and their effect on each of the
# ported Chmy-native computations. The load-bearing test is the NaN-containment one — that
# is the bug masking exists to fix, not a hypothetical.

@testset "ice masks" begin
    lx, ly, dx, dy = 8.0, 8.0, 1.0, 1.0

    # ice on the left half of the domain (x < 0), ice-free on the right
    left_half(x, y) = x < 0 ? 1000.0 : 0.0

    function setup(T = Float64; nlayers = nothing)
        grid = isnothing(nlayers) ? StaggeredGrid(T, lx, ly, dx, dy) :
               StaggeredGrid(T, lx, ly, dx, dy,
                             CorrectedVerticalLayering(T, QuadraticSigmaTransform(T, nlayers)))
        return grid, Runtime(grid), TopographicState(grid), MechanicState(grid)
    end

    @testset "icemasks!: is_ice and the neighbour ring" begin
        grid, rt, topo, _ = setup()
        fill_analytic!(topo.thickness.ice, rt.grid2d, (x, y) -> left_half(x, y))
        icemasks!(topo, rt)

        is_ice  = asarray(topo.mask.is_ice)[:, 2, 1]
        is_nbr  = asarray(topo.mask.is_ice_neighbour)[:, 2, 1]
        xc      = collect(grid.x)

        @test is_ice == (xc .< 0)
        # exactly one ice-free cell touches the ice, and it is the first one past x = 0
        @test count(is_nbr) == 1
        @test is_nbr[findfirst(>(0.0), xc)]
        # the ring is ice-*free* by construction: the two are disjoint
        @test !any(is_ice .& is_nbr)

        # is_ice_allowed is a user constraint, never derived — icemasks! leaves it alone
        @test all(iszero, asarray(topo.mask.is_ice_allowed))
    end

    @testset "icemasks!: H_min threshold" begin
        _, rt, topo, _ = setup()
        fill_analytic!(topo.thickness.ice, rt.grid2d, (x, y) -> x < 0 ? 1000.0 : 0.5)

        icemasks!(topo, rt)                  # H_min = 0: the 0.5 m film counts as ice
        @test all(asarray(topo.mask.is_ice))

        icemasks!(topo, rt; H_min = 1.0)     # ...and does not, above the threshold
        @test asarray(topo.mask.is_ice)[:, 2, 1] == (collect(rt.grid2d |> xcenters) .< 0)
    end

    # The node-class reduction: a face is active if *either* adjoining cell is, a corner if
    # any of four is. This encodes Chmy's convention that vertex `i` sits between centres
    # `i-1` and `i`, so it gets its own test rather than being trusted implicitly.
    @testset "node_active: the per-node-class reduction rules" begin
        grid, rt, topo, _ = setup()
        fill_analytic!(topo.thickness.ice, rt.grid2d, (x, y) -> left_half(x, y))
        icemasks!(topo, rt)
        mask = IceMask(topo.mask.is_ice)

        is_ice = asarray(topo.mask.is_ice)
        xc     = collect(grid.x)
        i_last = findlast(<(0.0), xc)          # last ice cell
        i_free = i_last + 1                    # first ice-free cell

        # aa: the cell itself
        @test node_active(mask, Pagos.NODE_AA, i_last, 2)
        @test !node_active(mask, Pagos.NODE_AA, i_free, 2)

        # acx: the face between cells i-1 and i. The face at i_free straddles the margin,
        # so it *is* active — this is what lets the margin advance.
        @test node_active(mask, Pagos.NODE_ACX, i_free, 2)
        @test is_ice[i_free - 1, 2, 1] && !is_ice[i_free, 2, 1]
        # ...while the next face out has ice on neither side
        @test !node_active(mask, Pagos.NODE_ACX, i_free + 1, 2)

        # acy: uniform in y here, so activity follows the cell's own column
        @test node_active(mask, Pagos.NODE_ACY, i_last, 2)
        @test !node_active(mask, Pagos.NODE_ACY, i_free, 2)

        # ab: any of the four surrounding cells
        @test node_active(mask, Pagos.NODE_AB, i_free, 2)
        @test !node_active(mask, Pagos.NODE_AB, i_free + 1, 2)

        # the z-Vertex variants share their horizontal activity
        @test node_active(mask, Pagos.NODE_ACX_AC, i_free, 2) ==
              node_active(mask, Pagos.NODE_ACX, i_free, 2)
        @test node_active(mask, Pagos.NODE_AA_AC, i_free, 2) ==
              node_active(mask, Pagos.NODE_AA, i_free, 2)
    end

    @testset "NoMask is active everywhere and compiles away" begin
        @test node_active(NoMask(), Pagos.NODE_AA, 1, 1)
        @test node_active(NoMask(), Pagos.NODE_AB, -5, 99)
        @test isbits(NoMask())
    end

    @testset "IceMask: the union of several fields widens the band" begin
        grid, rt, topo, _ = setup()
        fill_analytic!(topo.thickness.ice, rt.grid2d, (x, y) -> left_half(x, y))
        icemasks!(topo, rt)

        strict = IceMask(topo.mask.is_ice)
        wide   = IceMask(topo.mask.is_ice, topo.mask.is_ice_neighbour)
        i_free = findfirst(>(0.0), collect(grid.x))

        @test !node_active(strict, Pagos.NODE_AA, i_free, 2)
        @test node_active(wide, Pagos.NODE_AA, i_free, 2)
    end

    @testset "IceMask: `allowed` uses the stricter all-cells rule" begin
        grid, rt, topo, _ = setup()
        # ice everywhere, but only the left half is permitted
        fill_analytic!(topo.thickness.ice, rt.grid2d, (x, y) -> 1000.0)
        icemasks!(topo, rt)
        for k in -1:2, j in -1:(grid.ny + 2), i in -1:(grid.nx + 2)
            x, = coord(rt.grid2d, Pagos.NODE_AA, i, j, 1)
            topo.mask.is_ice_allowed[i, j, 1] = x < 0
        end
        mask = IceMask(topo.mask.is_ice; allowed = topo.mask.is_ice_allowed)

        xc     = collect(grid.x)
        i_last = findlast(<(0.0), xc)
        i_free = i_last + 1

        @test node_active(mask, Pagos.NODE_AA, i_last, 2)
        @test !node_active(mask, Pagos.NODE_AA, i_free, 2)
        # the face straddling the boundary has one disallowed side ⟹ blocked, unlike the
        # `active` rule which would keep it open
        @test !node_active(mask, Pagos.NODE_ACX, i_free, 2)
        @test node_active(IceMask(topo.mask.is_ice), Pagos.NODE_ACX, i_free, 2)
    end

    ###########################################################
    # The bug masking exists to fix
    ###########################################################

    # Verified before the mask existed: η is zero where there is no ice, `hlerp` averages
    # reciprocals, so τ_xy is NaN at ice-free faces — and the effective stress at `aa`
    # interpolates that back *one cell into the ice*, i.e. onto the margin.
    @testset "NaN containment: unmasked, NaN leaks one cell into the ice" begin
        grid, rt, _, mech = setup(; nlayers = 4)
        mat = MaterialState(grid)

        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> left_half(x, y))
        fill_analytic3d!(mat.eta_ice, rt.grid, (x, y, ζ) -> x < 0 ? 1.0e13 : 0.0)
        fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> 3.0 * y)
        fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> 0.0)
        fill_analytic3d!(mech.velocity.z, rt.grid, (x, y, ζ) -> 0.0)

        mb = BlatterPattynMomentumBalance()
        xc = collect(grid.x)
        i_last = findlast(<(0.0), xc)

        # --- unmasked: the leak ---
        raw_strainrate!(mech, mb, rt)
        deviatoric_stress!(mech, mat, rt)
        eff = asarray(mech.stress.effective)
        @test Pagos.hasnan(mech.stress.xy)
        @test any(isnan, eff[i_last, :, :])          # NaN on an ICE-COVERED cell

        # --- masked: contained. A fresh state, so a leftover NaN from the run above
        # cannot be what makes this pass or fail. ---
        _, rt2, topo, mech2 = setup(; nlayers = 4)
        mat2 = MaterialState(grid)
        fill_analytic!(topo.thickness.ice, rt2.grid2d, (x, y) -> left_half(x, y))
        fill_analytic!(mech2.topography.thickness, rt2.grid2d, (x, y) -> left_half(x, y))
        fill_analytic3d!(mat2.eta_ice, rt2.grid, (x, y, ζ) -> x < 0 ? 1.0e13 : 0.0)
        fill_analytic3d!(mech2.velocity.x, rt2.grid, (x, y, ζ) -> 3.0 * y)
        fill_analytic3d!(mech2.velocity.y, rt2.grid, (x, y, ζ) -> 0.0)
        fill_analytic3d!(mech2.velocity.z, rt2.grid, (x, y, ζ) -> 0.0)
        icemasks!(topo, rt2)
        mask = IceMask(topo.mask.is_ice)

        raw_strainrate!(mech2, mb, rt2, mask)
        deviatoric_stress!(mech2, mat2, rt2, mask)
        mech = mech2

        @test !Pagos.hasnan(mech.stress.xy)
        @test !Pagos.hasnan(mech.stress.effective)
        @test all(isfinite, asarray(mech.stress.effective))
        # ice-covered cells still carry a stress; ice-free ones are exactly zero
        @test asarray(mech.stress.effective)[i_last, 2, 1] > 0
        @test all(iszero, asarray(mech.stress.effective)[(i_last + 1):end, :, :])
    end

    # `txy`/`txz`/`tyz` are three separate hardcoded call sites in the kernel, each naming
    # its own node class (`NODE_AB`/`NODE_ACX_AC`/`NODE_ACY_AC`). The test above only ever
    # gives `xy` a nonzero, NaN-risking value (velocity varies in y only, not with ζ), so a
    # copy-paste slip in the `xz`/`yz` lines — e.g. gating `txz` on `NODE_AB` instead of
    # `NODE_ACX_AC` — would go uncaught: `_mask_cells` reduces those two classes over the
    # same *horizontal* footprint as `acx`, so the mistake would still compile and mostly
    # look right, just at the wrong vertical interface. This variant gives the velocity
    # vertical shear instead, so `ε̇_xz`/`τ_xz` are the nonzero, NaN-risking components.
    @testset "NaN containment: the vertical-shear branch independently" begin
        grid, rt, _, mech = setup(; nlayers = 4)
        mat = MaterialState(grid)

        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> left_half(x, y))
        fill_analytic3d!(mat.eta_ice, rt.grid, (x, y, ζ) -> x < 0 ? 1.0e13 : 0.0)
        fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> 5.0 * ζ)   # ∂u/∂z ≠ 0
        fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> 0.0)
        fill_analytic3d!(mech.velocity.z, rt.grid, (x, y, ζ) -> 0.0)

        mb = BlatterPattynMomentumBalance()

        # --- unmasked: xz (not xy) is where the NaN shows up here ---
        raw_strainrate!(mech, mb, rt)
        @test !Pagos.hasnan(mech.strainrate.xz)      # the strain rate itself is finite
        deviatoric_stress!(mech, mat, rt)
        @test Pagos.hasnan(mech.stress.xz)           # NaN * 0 == NaN, in the xz branch

        # --- masked: contained, and the ice-covered value survives ---
        _, rt2, topo, mech2 = setup(; nlayers = 4)
        mat2 = MaterialState(grid)
        fill_analytic!(topo.thickness.ice, rt2.grid2d, (x, y) -> left_half(x, y))
        fill_analytic!(mech2.topography.thickness, rt2.grid2d, (x, y) -> left_half(x, y))
        fill_analytic3d!(mat2.eta_ice, rt2.grid, (x, y, ζ) -> x < 0 ? 1.0e13 : 0.0)
        fill_analytic3d!(mech2.velocity.x, rt2.grid, (x, y, ζ) -> 5.0 * ζ)
        fill_analytic3d!(mech2.velocity.y, rt2.grid, (x, y, ζ) -> 0.0)
        fill_analytic3d!(mech2.velocity.z, rt2.grid, (x, y, ζ) -> 0.0)
        icemasks!(topo, rt2)
        mask = IceMask(topo.mask.is_ice)

        raw_strainrate!(mech2, mb, rt2, mask)
        deviatoric_stress!(mech2, mat2, rt2, mask)

        xc     = collect(grid.x)
        i_last = findlast(<(0.0), xc)
        @test !Pagos.hasnan(mech2.stress.xz)
        @test interior(mech2.stress.xz)[i_last, 2, 3] != 0.0      # ice-covered: real value
        @test all(iszero, interior(mech2.stress.xz)[(i_last + 2):end, :, :])  # clear of ice
    end

    # Masking the strain rate alone is *not* enough, because NaN * 0 == NaN. This pins the
    # reason the mask has to be threaded into the stress kernel too.
    @testset "masking the strain rate alone does not contain the NaN" begin
        grid, rt, topo, mech = setup(; nlayers = 4)
        mat = MaterialState(grid)

        fill_analytic3d!(mat.eta_ice, rt.grid, (x, y, ζ) -> x < 0 ? 1.0e13 : 0.0)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> left_half(x, y))
        fill_analytic!(topo.thickness.ice, rt.grid2d, (x, y) -> left_half(x, y))
        fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> 0.0)
        fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> 0.0)
        fill_analytic3d!(mech.velocity.z, rt.grid, (x, y, ζ) -> 0.0)
        icemasks!(topo, rt)
        mask = IceMask(topo.mask.is_ice)

        # strain rate masked (so ε̇_xy is exactly 0 on ice-free nodes) but stress NOT masked
        raw_strainrate!(mech, BlatterPattynMomentumBalance(), rt, mask)
        @test !Pagos.hasnan(mech.strainrate.xy)
        @test all(iszero, asarray(mech.strainrate.xy))

        deviatoric_stress!(mech, mat, rt)            # unmasked: NaN * 0 == NaN
        @test Pagos.hasnan(mech.stress.xy)
    end

    ###########################################################
    # drivingstress!/surface_gradient!/velocitygradients!: masked directly
    ###########################################################
    #
    # These three all use the permissive rule (no hlerp involved, so no NaN risk), and
    # their masking was previously only exercised *transitively*, through raw_strainrate!'s
    # composition in the NaN-containment test above, or with NoMask() in the equivalence
    # test below — never with a real IceMask checked against its own output. Manually
    # verified against the code before writing these (a first draft of this check had a
    # bug of its own: it filled `mech.topography.surface` while calling
    # `surface_gradient!(topo, ...)`, which reads `topo.elevation.surface` — a field that
    # was never filled, so every result was a tautological zero. Caught by checking the
    # *inactive* value matched *and differed from* the active one, not just its sign.)

    @testset "surface_gradient! and drivingstress!: masked directly" begin
        grid, rt, topo, mech = setup()
        cst = Constants{Float64}()

        fill_analytic!(topo.elevation.surface, rt.grid2d, (x, y) -> 0.02x)
        fill_analytic!(mech.topography.surface, rt.grid2d, (x, y) -> 0.02x)
        fill_analytic!(topo.thickness.ice, rt.grid2d, (x, y) -> left_half(x, y))
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> left_half(x, y))
        icemasks!(topo, rt)
        mask = IceMask(topo.mask.is_ice)

        surface_gradient!(topo, rt, mask)
        drivingstress!(mech, cst, rt, mask)

        xc     = collect(grid.x)
        i_free = findfirst(>(0.0), xc)     # the face straddling the margin
        i_far  = i_free + 1                # a face with ice on neither side

        # Both fields at the margin face: real and nonzero, matching what the (unmasked)
        # formula gives there — the permissive rule keeps this face open, so masking must
        # not perturb its result. `H` at the face is `lerp`-averaged across the margin
        # (`left_half` gives 1000 on one side, 0 on the other ⟹ 500 at the face), not the
        # full interior thickness.
        @test interior(topo.elevation.surface_dx)[i_free, 2, 1] ≈ 0.02
        @test interior(mech.stress.driving_x)[i_free, 2, 1] ≈
              cst.density_ice * cst.gravity * 500.0 * 0.02

        # ...while a face with no ice on either side is exactly zero, not just small.
        @test interior(topo.elevation.surface_dx)[i_far, 2, 1] == 0.0
        @test interior(mech.stress.driving_x)[i_far, 2, 1] == 0.0
        @test interior(topo.elevation.surface_dx)[i_free, 2, 1] != 0.0   # the two differ
    end

    # `velocitygradients!` is exercised transitively above (through `raw_strainrate!`'s
    # state-level composition), but never checked on its own output directly with a real
    # mask. This pins that the gradient fields themselves — not just the downstream strain
    # rate — are exactly zero off the active region.
    @testset "velocitygradients! masked directly" begin
        grid, rt, topo, mech = setup(; nlayers = 4)
        fill_analytic!(topo.thickness.ice, rt.grid2d, (x, y) -> left_half(x, y))
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> left_half(x, y))
        fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> 3.0 * y)
        fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> 0.0)
        fill_analytic3d!(mech.velocity.z, rt.grid, (x, y, ζ) -> 0.0)
        icemasks!(topo, rt)
        mask = IceMask(topo.mask.is_ice)

        velocitygradients!(mech.velocity, mech.topography.thickness, rt, mask)

        xc     = collect(grid.x)
        i_last = findlast(<(0.0), xc)
        i_far  = i_last + 2       # two cells past the margin: no adjoining ice anywhere

        @test interior(mech.velocity.x_dy)[i_last, 2, 1] ≈ 3.0     # active: real value
        @test interior(mech.velocity.x_dy)[i_far, 2, 1] == 0.0     # inactive: exact zero
    end

    ###########################################################
    # Advection: conservation and margin advance
    ###########################################################

    # The reason the mask goes on the fluxes and not the tendency: with flux masking the
    # divergence still telescopes, so total mass change equals the net boundary flux plus
    # the mass balance, to machine precision.
    @testset "flux masking preserves exact conservation" begin
        grid, rt, topo, mech = setup()
        fill_analytic!(topo.thickness.ice, rt.grid2d,
                       (x, y) -> x < 0 ? 500 + 100sin(3x) * cos(2y) : 0.0)
        fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d, (x, y) -> sin(x) + 0.5cos(y))
        fill_analytic!(mech.velocity.depthaverage_y, rt.grid2d, (x, y) -> cos(2x) - 0.3sin(y))
        icemasks!(topo, rt)
        mask = IceMask(topo.mask.is_ice)

        advect!(topo, mech, CenteredAdvection(), rt, mask)

        qx, qy = asarray(mech.flux.x), asarray(mech.flux.y)
        influx = (sum(qx[1, :, :]) - sum(qx[end, :, :])) * dy +
                 (sum(qy[:, 1, :]) - sum(qy[:, end, :])) * dx
        @test sum(asarray(topo.thickness.ice_dt)) * dx * dy ≈ influx rtol = 1e-12
    end

    # The margin must still be able to advance: the face between the last ice cell and the
    # first empty one carries flux, so the empty cell gets a positive tendency.
    @testset "the margin can still advance under an is_ice mask" begin
        grid, rt, topo, mech = setup()
        fill_analytic!(topo.thickness.ice, rt.grid2d, (x, y) -> left_half(x, y))
        fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d, (x, y) -> 1.0)   # eastward
        fill_analytic!(mech.velocity.depthaverage_y, rt.grid2d, (x, y) -> 0.0)
        icemasks!(topo, rt)

        advect!(topo, mech, UpwindAdvection(), rt, IceMask(topo.mask.is_ice))

        xc     = collect(grid.x)
        i_free = findfirst(>(0.0), xc)
        @test asarray(mech.flux.x)[i_free, 2, 1] > 0        # flux into the empty cell
        @test asarray(topo.thickness.ice_dt)[i_free, 2, 1] > 0   # ...so it gains ice
        # ...and nothing is happening two cells out
        @test asarray(mech.flux.x)[i_free + 1, 2, 1] == 0
    end

    # A hard `allowed` constraint blocks advance instead, and does so conservatively: the
    # ice piles up in the donor cell rather than vanishing.
    @testset "an `allowed` wall blocks advance without destroying mass" begin
        grid, rt, topo, mech = setup()
        fill_analytic!(topo.thickness.ice, rt.grid2d, (x, y) -> left_half(x, y))
        fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d, (x, y) -> 1.0)
        fill_analytic!(mech.velocity.depthaverage_y, rt.grid2d, (x, y) -> 0.0)
        icemasks!(topo, rt)
        for k in -1:2, j in -1:(grid.ny + 2), i in -1:(grid.nx + 2)
            x, = coord(rt.grid2d, Pagos.NODE_AA, i, j, 1)
            topo.mask.is_ice_allowed[i, j, 1] = x < 0
        end
        mask = IceMask(topo.mask.is_ice; allowed = topo.mask.is_ice_allowed)

        advect!(topo, mech, UpwindAdvection(), rt, mask)

        xc     = collect(grid.x)
        i_last = findlast(<(0.0), xc)
        @test asarray(mech.flux.x)[i_last + 1, 2, 1] == 0        # the wall
        @test asarray(topo.thickness.ice_dt)[i_last + 1, 2, 1] == 0  # no advance
        @test asarray(topo.thickness.ice_dt)[i_last, 2, 1] > 0   # ice piles up instead

        # ...and mass is still conserved: what the wall blocks stays in the domain
        qx, qy = asarray(mech.flux.x), asarray(mech.flux.y)
        influx = (sum(qx[1, :, :]) - sum(qx[end, :, :])) * dy +
                 (sum(qy[:, 1, :]) - sum(qy[:, end, :])) * dx
        @test sum(asarray(topo.thickness.ice_dt)) * dx * dy ≈ influx rtol = 1e-12
    end

    ###########################################################
    # The unmasked path must be unchanged
    ###########################################################

    # `NoMask()` is the default and dispatches at compile time, so every existing result has
    # to be bit-for-bit what it was before mask support existed.
    @testset "NoMask is bit-for-bit identical to passing no mask at all" begin
        cst = Constants{Float64}()

        function run(args...)
            grid, rt, topo, mech = setup(; nlayers = 4)
            mat = MaterialState(grid)
            fill_analytic!(mech.topography.surface, rt.grid2d, (x, y) -> 0.01x^2 - 0.02y)
            fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> 800 + 10x)
            fill_analytic!(topo.thickness.ice, rt.grid2d, (x, y) -> 800 + 10x)
            fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> 2x + 3y)
            fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> -x + 0.5y)
            fill_analytic3d!(mech.velocity.z, rt.grid, (x, y, ζ) -> 0.0)
            fill_analytic3d!(mat.eta_ice, rt.grid, (x, y, ζ) -> 1.0e13)
            fill_analytic!(mech.velocity.depthaverage_x, rt.grid2d, (x, y) -> 2.0)
            fill_analytic!(mech.velocity.depthaverage_y, rt.grid2d, (x, y) -> 1.0)

            drivingstress!(mech, cst, rt, args...)
            surface_gradient!(topo, rt, args...)
            raw_strainrate!(mech, BlatterPattynMomentumBalance(), rt, args...)
            deviatoric_stress!(mech, mat, rt, args...)
            advect!(topo, mech, CenteredAdvection(), rt, args...)
            return (copy(asarray(mech.stress.driving_x)),
                    copy(asarray(topo.elevation.surface_dx)),
                    copy(asarray(mech.strainrate.xy)),
                    copy(asarray(mech.strainrate.effective)),
                    copy(asarray(mech.stress.xy)),
                    copy(asarray(mech.stress.effective)),
                    copy(asarray(mech.flux.x)),
                    copy(asarray(topo.thickness.ice_dt)))
        end

        default  = run()
        explicit = run(NoMask())
        for (a, b) in zip(default, explicit)
            @test a == b                # bit-for-bit, not ≈
        end
    end

    # -------------------------------------------------------------------
    # momentum_mask!: connectivity to grounded ice. The reference is a serial flood fill —
    # the same oracle the AIS geometry was validated against — because the kernel's
    # iterative propagation is the thing under test, not the definition.
    # -------------------------------------------------------------------

    # Deliberately not `momentum_mask!`'s algorithm: an explicit stack-based flood fill,
    # so agreement is evidence about the propagation kernel rather than a tautology.
    function floodfill(ice, grounded)
        nx, ny = size(ice)
        out = ice .& grounded
        stack = [(i, j) for j in 1:ny, i in 1:nx if out[i, j]]
        while !isempty(stack)
            i, j = pop!(stack)
            for (a, b) in ((i-1, j), (i+1, j), (i, j-1), (i, j+1))
                if 1 <= a <= nx && 1 <= b <= ny && ice[a, b] && !out[a, b]
                    out[a, b] = true
                    push!(stack, (a, b))
                end
            end
        end
        return out
    end

    function momentum_setup(H, grounded)
        nx, ny = size(H)
        g  = StaggeredGrid(Float64, nx * 1.0, ny * 1.0, 1.0, 1.0)
        r  = Runtime(g)
        tp = TopographicState(g)
        setdata!(tp.thickness.ice, H)
        setdata!(tp.mask.is_grounded, grounded)
        icemasks!(tp, r)
        momentum_mask!(tp, r)
        return asarray(tp.mask.is_momentum_solved)[:, :, 1]
    end

    @testset "momentum_mask!: excludes detached ice, keeps everything connected" begin
        nx, ny = 12, 12
        H = zeros(nx, ny)
        H[1:5, :]   .= 1000.0        # grounded sheet + attached shelf
        H[9:10, 3:4] .= 300.0        # detached berg
        H[7, 7]      = 200.0         # single detached cell
        grounded = falses(nx, ny)
        grounded[1:2, :] .= true

        solved = momentum_setup(H, grounded)

        @test solved == floodfill(H .> 0, grounded)
        @test all(solved[1:5, :])                     # sheet and its shelf: solved
        @test !any(solved[9:10, 3:4])                 # berg: not solved
        @test !solved[7, 7]
        @test !any(solved[H .== 0])                   # never marks ice-free cells

        # A berg that touches the sheet only at a corner stays excluded: our membrane
        # stress gates the `ab` corner term on `node_fully_active`, so a diagonal contact
        # transmits exactly zero stress. 4-connectivity is the discretization's rule, not
        # a convention — 8-connectivity would mark this cell load-bearing.
        H_diag = zeros(nx, ny)
        H_diag[1:5, 1:5] .= 1000.0
        H_diag[6, 6] = 300.0                          # touches (5,5) only diagonally
        solved_diag = momentum_setup(H_diag, grounded)
        @test all(solved_diag[1:5, 1:5])
        @test !solved_diag[6, 6]
    end

    # The seed is `is_ice & is_grounded`, so grounded-but-ice-free bedrock must not seed a
    # component, and floating ice reachable only through such a cell stays excluded.
    @testset "momentum_mask!: grounded but ice-free cells do not seed" begin
        nx, ny = 8, 8
        H = zeros(nx, ny)
        H[5:6, 5:6] .= 400.0                          # floating patch, no grounded ice at all
        grounded = falses(nx, ny)
        grounded[1:2, 1:2] .= true                    # grounded, but H = 0 there

        solved = momentum_setup(H, grounded)
        @test !any(solved)                            # nothing is well-posed: no seed
        @test solved == floodfill(H .> 0, grounded)
    end

    # A serpentine channel is the case where propagation needs many sweeps: the path length
    # is far longer than the domain diameter, which is exactly what `maxsweeps` guards.
    @testset "momentum_mask!: serpentine path, and maxsweeps throws rather than truncating" begin
        nx, ny = 21, 21
        H = zeros(nx, ny)
        for (r, row) in enumerate(1:2:ny)             # horizontal bars ...
            H[:, row] .= 500.0
            r <= (ny - 1) ÷ 2 &&                      # ... joined alternately at each end
                (H[isodd(r) ? nx : 1, row + 1] = 500.0)
        end
        grounded = falses(nx, ny)
        grounded[1, 1] = true

        solved = momentum_setup(H, grounded)
        @test solved == floodfill(H .> 0, grounded)
        @test all(solved[:, 1])
        @test all(solved[:, ny])                      # reached the far end of the serpentine

        # Truncating would silently freeze connected ice at zero velocity, so it throws.
        g  = StaggeredGrid(Float64, nx * 1.0, ny * 1.0, 1.0, 1.0)
        r  = Runtime(g)
        tp = TopographicState(g)
        setdata!(tp.thickness.ice, H)
        setdata!(tp.mask.is_grounded, grounded)
        icemasks!(tp, r)
        @test_throws ErrorException momentum_mask!(tp, r; maxsweeps = 1)
    end

    # The point of the whole exercise: swap the mask, no solver change.
    @testset "momentum_mask! feeds IceMask directly" begin
        nx, ny = 10, 10
        H = zeros(nx, ny); H[1:4, :] .= 1000.0; H[8:9, 4:5] .= 200.0
        grounded = falses(nx, ny); grounded[1:2, :] .= true

        g  = StaggeredGrid(Float64, nx * 1.0, ny * 1.0, 1.0, 1.0)
        r  = Runtime(g)
        tp = TopographicState(g)
        setdata!(tp.thickness.ice, H)
        setdata!(tp.mask.is_grounded, grounded)
        icemasks!(tp, r)
        momentum_mask!(tp, r)

        mask = IceMask(tp.mask.is_momentum_solved)
        @test node_active(mask, Pagos.NODE_AA, 2, 5)          # inside the sheet
        @test !node_active(mask, Pagos.NODE_AA, 9, 5)         # inside the berg
        # A face is permissive (any adjacent cell active), the `ab` corner strict — the
        # existing IceMask rules apply unchanged to this field, which is the design point.
        @test node_active(mask, Pagos.NODE_ACX, 5, 5)         # sheet edge face: still active
        @test !node_fully_active(mask, Pagos.NODE_ACX, 5, 5)
    end
end

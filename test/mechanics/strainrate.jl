using Pagos
using Test

# nz > 1 grid (the public constructor fixes nz = 1); tensor fields are always 3D.
function grid3d(T, lx, ly, dx, dy, nz)
    g = RegularGrid(T, lx, ly, dx, dy)
    z = collect(range(T(0), T(1); length = nz))
    return RegularGrid(g.nx, g.ny, nz, g.x, g.y, z, g.dx, g.dy, g.dz,
        g.Lon, g.Lat, g.area, g.distortion, g.basins, g.regions)
end

# `raw_strainrate!(strainrate, velocity, momentum)` is the low-level entry point that
# reconstructs/dispatches ε̇_zz exactly like the old `strainrate!(mech[, balance])` used
# to (see raw_strainrate! docstring). `SSAMomentumBalance` stands in for "any
# depth-integrated balance" below: the fallback method is shared by all of them.

# ---------------------------------------------------------------------------
# Plane dynamics (nz == 1): vertical velocity gradients are zero
# ---------------------------------------------------------------------------

@testset "strainrate! — plane (nz = 1)" begin
    grid = RegularGrid(Float64, 4.0, 4.0, 1.0, 1.0)
    mech = MechanicState(grid)

    fill!(mech.velocity.x_dx,  2.0)   # ε̇xx = 2
    fill!(mech.velocity.y_dy, -1.0)   # ε̇yy = -1
    fill!(mech.velocity.x_dy,  1.0)
    fill!(mech.velocity.y_dx,  3.0)   # ε̇xy = (1 + 3)/2 = 2

    raw_strainrate!(mech.strainrate, mech.velocity, SSAMomentumBalance())

    @test all(mech.strainrate.xx .≈  2.0)
    @test all(mech.strainrate.yy .≈ -1.0)
    @test all(mech.strainrate.xy .≈  2.0)
    # ε̇zz reconstructed: -(2 - 1) = -1
    @test all(mech.strainrate.zz .≈ -1.0)
    # no resolved vertical shear
    @test all(mech.strainrate.xz .== 0.0)
    @test all(mech.strainrate.yz .== 0.0)
    # effective: √(½(4 + 1 + 1) + 4) = √(3 + 4) = √7
    @test all(mech.strainrate.effective .≈ sqrt(7.0))
end

# ---------------------------------------------------------------------------
# Full-column dynamics (nz > 1): vertical shear present
# ---------------------------------------------------------------------------

@testset "strainrate! — full-column (nz > 1)" begin
    grid = grid3d(Float64, 4.0, 4.0, 1.0, 1.0, 4)
    mech = MechanicState(grid)

    fill!(mech.velocity.x_dx,  2.0)   # ε̇xx = 2
    fill!(mech.velocity.y_dy, -1.0)   # ε̇yy = -1
    fill!(mech.velocity.x_dy,  1.0)
    fill!(mech.velocity.y_dx,  3.0)   # ε̇xy = 2
    fill!(mech.velocity.x_dz,  4.0)
    fill!(mech.velocity.z_dx,  0.0)   # ε̇xz = (4 + 0)/2 = 2
    fill!(mech.velocity.y_dz, -2.0)
    fill!(mech.velocity.z_dy,  2.0)   # ε̇yz = (-2 + 2)/2 = 0

    raw_strainrate!(mech.strainrate, mech.velocity, SSAMomentumBalance())

    @test all(mech.strainrate.xx .≈  2.0)
    @test all(mech.strainrate.yy .≈ -1.0)
    @test all(mech.strainrate.zz .≈ -1.0)   # reconstructed
    @test all(mech.strainrate.xy .≈  2.0)
    @test all(mech.strainrate.xz .≈  2.0)
    @test all(mech.strainrate.yz .≈  0.0)
    # effective: √(½(4 + 1 + 1) + 4 + 4 + 0) = √(3 + 8) = √11
    @test all(mech.strainrate.effective .≈ sqrt(11.0))
end

# ---------------------------------------------------------------------------
# MomentumBalance dispatch: full-column uses ∂w/∂z, others reconstruct ε̇zz
# ---------------------------------------------------------------------------

@testset "strainrate! — MomentumBalance dispatch on ε̇zz" begin
    grid = grid3d(Float64, 4.0, 4.0, 1.0, 1.0, 4)
    mech = MechanicState(grid)

    fill!(mech.velocity.x_dx,  2.0)   # ε̇xx = 2
    fill!(mech.velocity.y_dy, -1.0)   # ε̇yy = -1
    # ∂w/∂z deliberately ≠ -(ε̇xx + ε̇yy) = -1, to distinguish the two paths
    fill!(mech.velocity.z_dz,  5.0)

    # Full-column balances read ∂w/∂z directly
    for balance in (StokesMomentumBalance(), BlatterPattynMomentumBalance())
        fill!(mech.strainrate.zz, 0.0)
        raw_strainrate!(mech.strainrate, mech.velocity, balance)
        @test all(mech.strainrate.zz .≈ 5.0)
        # only the diagonal is set: √(½(4 + 1 + 25)) = √15
        @test all(mech.strainrate.effective .≈ sqrt(15.0))
    end

    # Depth-integrated balance reconstructs from incompressibility
    raw_strainrate!(mech.strainrate, mech.velocity, SSAMomentumBalance())
    @test all(mech.strainrate.zz .≈ -1.0)
    # ... any other depth-integrated balance hits the same generic fallback
    raw_strainrate!(mech.strainrate, mech.velocity, SIAMomentumBalance())
    @test all(mech.strainrate.zz .≈ -1.0)
end

# ---------------------------------------------------------------------------
# Strain rate → stress consistency: τ_ij = 2 η ε̇_ij through the two kernels
# ---------------------------------------------------------------------------

@testset "strainrate! feeds deviatoric_stress!" begin
    grid = grid3d(Float64, 5.0, 5.0, 1.0, 1.0, 3)
    mech = MechanicState(grid)
    mat  = MaterialState(grid)

    fill!(mech.velocity.x_dx,  2.0); fill!(mech.velocity.y_dy, -1.0)
    fill!(mech.velocity.x_dy,  1.0); fill!(mech.velocity.y_dx,  3.0)
    fill!(mech.velocity.x_dz,  4.0); fill!(mech.velocity.z_dx,  0.0)
    fill!(mech.velocity.y_dz, -2.0); fill!(mech.velocity.z_dy,  2.0)
    η = 2.0
    fill!(mat.eta_ice, η)

    raw_strainrate!(mech.strainrate, mech.velocity, SSAMomentumBalance())
    deviatoric_stress!(mech, mat)

    @test all(mech.stress.xx .≈ 2η .* mech.strainrate.xx)
    @test all(mech.stress.xz .≈ 2η .* mech.strainrate.xz)
    # effective stress is 2 η times effective strain rate
    @test all(mech.stress.effective .≈ 2η .* mech.strainrate.effective)
end

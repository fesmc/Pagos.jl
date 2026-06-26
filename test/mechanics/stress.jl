using Pagos
using Test

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# The public RegularGrid constructor fixes nz = 1. For the full-column case we
# reconstruct it with nz > 1 and a matching z vector. The tensor/stress fields are
# always 3D (nx, ny, nz); depth-averaged dynamics simply use nz == 1.
function grid3d(T, lx, ly, dx, dy, nz)
    g = RegularGrid(T, lx, ly, dx, dy)
    z = collect(range(T(0), T(1); length = nz))
    return RegularGrid(g.nx, g.ny, nz, g.x, g.y, z, g.dx, g.dy, g.dz,
        g.Lon, g.Lat, g.area, g.distortion, g.basins, g.regions)
end

# ---------------------------------------------------------------------------
# Depth-averaged dynamics (nz == 1): vertical-shear strain rates are zero
# ---------------------------------------------------------------------------

@testset "deviatoric_stress! — depth-averaged (nz = 1)" begin
    grid = RegularGrid(Float64, 4.0, 4.0, 1.0, 1.0)
    mech = MechanicState(grid)
    mat  = MaterialState(grid)

    @test ndims(mech.strain_rate_dxx) == 3   # tensor fields are always 3D
    @test size(mech.strain_rate_dxx, 3) == 1 # one layer for depth-averaged dynamics

    # Constant viscosity and constant (uniform) strain rate. The vertical-shear
    # strain rates are left at zero (no resolved vertical shear in 2D dynamics).
    η = 2.0
    fill!(mat.eta_ice, η)
    fill!(mech.strain_rate_dxx,  3.0)
    fill!(mech.strain_rate_dyy, -1.0)
    fill!(mech.strain_rate_dxy,  0.5)

    deviatoric_stress!(mech, mat)

    # τ_ij = 2 η ε̇_ij
    @test all(mech.stress_xx .≈ 2 * η *  3.0)   # 12
    @test all(mech.stress_yy .≈ 2 * η * -1.0)   # -4
    @test all(mech.stress_xy .≈ 2 * η *  0.5)   # 2

    # τ_zz reconstructed from the traceless identity: -(τxx + τyy) = -(12 - 4) = -8.
    @test all(mech.stress_zz .≈ -(12.0 - 4.0))  # -8
    # Vertical-shear stresses vanish in the depth-averaged case.
    @test all(mech.stress_xz .== 0.0)
    @test all(mech.stress_yz .== 0.0)

    # Effective stress (unified formula): √(½(144 + 16 + 64) + 4) = √(112 + 4) = √116.
    # Equivalently √(τxx² + τyy² + τxx τyy + τxy²) with the reconstructed τzz.
    @test all(mech.stress_effective .≈ sqrt(116.0))
end

# ---------------------------------------------------------------------------
# Full-column dynamics (nz > 1)
# ---------------------------------------------------------------------------

@testset "deviatoric_stress! — full-column (nz > 1)" begin
    grid = grid3d(Float64, 4.0, 4.0, 1.0, 1.0, 4)
    mech = MechanicState(grid)
    mat  = MaterialState(grid)

    @test size(mech.strain_rate_dxx, 3) == 4

    # Incompressible strain rate (ε̇xx + ε̇yy + ε̇zz = 3 - 1 - 2 = 0), so the
    # reconstructed τ_zz = -(τxx + τyy) matches the direct 2 η ε̇_zz value.
    η = 2.0
    fill!(mat.eta_ice, η)
    fill!(mech.strain_rate_dxx,  3.0)
    fill!(mech.strain_rate_dyy, -1.0)
    fill!(mech.strain_rate_dzz, -2.0)
    fill!(mech.strain_rate_dxy,  0.5)
    fill!(mech.strain_rate_dxz,  1.0)
    fill!(mech.strain_rate_dyz, -0.5)

    deviatoric_stress!(mech, mat)

    @test all(mech.stress_xx .≈ 2 * η *  3.0)   # 12
    @test all(mech.stress_yy .≈ 2 * η * -1.0)   # -4
    @test all(mech.stress_zz .≈ -(12.0 - 4.0))  # -8, reconstructed = 2 η ε̇_zz here
    @test all(mech.stress_xy .≈ 2 * η *  0.5)   # 2
    @test all(mech.stress_xz .≈ 2 * η *  1.0)   # 4
    @test all(mech.stress_yz .≈ 2 * η * -0.5)   # -2

    # Effective stress: √(½(144 + 16 + 64) + 4 + 16 + 4) = √(112 + 24) = √136
    @test all(mech.stress_effective .≈ sqrt(136.0))
end

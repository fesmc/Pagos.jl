# ---------------------------------------------------------------------------
# Stress terms (`src/mechanics/stress.jl`).
#
# The other three terms of the force balance: the driving stress (computed once, before the
# PT loop), the basal stress (recomputed every iteration under `ActiveFrictionUpdate`), and
# the full 3D deviatoric tensor (the Blatter–Pattyn/Stokes path, and what a thermodynamic
# coupling reads for strain heating).
# ---------------------------------------------------------------------------

"""
    stress_suite(fx)

Benchmarks for the driving stress and its surface-gradient half, the basal stress, and the
fused 3D deviatoric-stress kernel.
"""
function stress_suite(fx)
    (; rt, mech, mat, cst, mask, backend) = fx
    g = BenchmarkGroup()

    # State-level `drivingstress!` is `surface_gradient!` plus the ρgH scaling; both are
    # tracked so a regression lands on the stencil or on the scaling, not on their sum.
    g["drivingstress!"] = @benchmarkable begin
        drivingstress!($mech, $cst, $rt, $mask)
        sync!($backend)
    end
    g["surface_gradient!"] = @benchmarkable begin
        surface_gradient!($mech.stress.driving_x, $mech.stress.driving_y,
                          $mech.topography.surface, $rt, $mask)
        sync!($backend)
    end

    # Per PT iteration under the default `ActiveFrictionUpdate` — cheap per call, but paid
    # as many times as the solve takes iterations.
    g["basalstress!"] = @benchmarkable begin
        basalstress!($mech, $rt, $mask)
        sync!($backend)
    end

    # Two launches over the full column (the tensor, then its second invariant). Not on the
    # SSA/DIVA path — this is what the higher-order balances and the strain-heating
    # coupling pay.
    g["deviatoric_stress!"] = @benchmarkable begin
        deviatoric_stress!($mech, $mat, $rt, $mask)
        sync!($backend)
    end

    return g
end

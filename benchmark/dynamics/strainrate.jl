# ---------------------------------------------------------------------------
# Velocity gradients, strain rates and the membrane stress
# (`src/mechanics/strainrate.jl`).
#
# The first half of every PT iteration: `pseudo_rate!` calls
# `depthaverage_velocitygradients!` → `effective_strainrate_*!` → `membranestress!` before
# anything else happens. Timed individually here so a regression can be attributed to one
# kernel instead of to "the solve got slower".
#
# Each of these is a pure function of state the fixture already holds, and each writes into
# fields nothing else in this group reads, so they need no `setup` — repeated evaluations
# recompute the same answer from the same inputs.
# ---------------------------------------------------------------------------

"""
    strainrate_suite(fx)

Benchmarks for the Chmy-native (C-grid staggered) velocity-gradient and strain-rate
kernels, plus the SSA/DIVA membrane stress they feed.
"""
function strainrate_suite(fx)
    (; rt, mech, mask, backend) = fx
    g = BenchmarkGroup()

    ssa  = SSAMomentumBalance()
    diva = DIVAMomentumBalance()

    # Full-column gradients (3D launch) vs. the depth-averaged four (2D launch): the PT
    # loop only pays the latter, so the gap between them is the cost DIVA's 3D
    # reconstruction adds when `velocities3D!` is finally called.
    g["velocitygradients!"] = @benchmarkable begin
        velocitygradients!($mech.velocity, $mech.topography.thickness, $rt, $mask)
        sync!($backend)
    end
    g["depthaverage_velocitygradients!"] = @benchmarkable begin
        depthaverage_velocitygradients!($mech.velocity, $rt, $mask)
        sync!($backend)
    end

    # State-level `raw_strainrate!` chains three launches (gradients → tensor → invariant);
    # tracked as the one call a consumer makes, since that is the unit anything outside
    # `strainrate.jl` uses.
    g["raw_strainrate! (state)"] = @benchmarkable begin
        raw_strainrate!($mech, $diva, $rt, $mask)
        sync!($backend)
    end

    # SSA's invariant is depth-independent (2D launch, Eq. 12); DIVA's adds the
    # vertical-shear terms and runs per layer (3D launch, Eq. 13). Both are on the PT hot
    # path for their respective balance, so both are tracked.
    g["effective_strainrate_ssa!"] = @benchmarkable begin
        effective_strainrate_ssa!($mech.strainrate, $mech.velocity, $rt, $mask)
        sync!($backend)
    end
    g["effective_strainrate_diva!"] = @benchmarkable begin
        effective_strainrate_diva!($mech, $rt, $mask)
        sync!($backend)
    end

    g["membranestress! (SSA)"] = @benchmarkable begin
        membranestress!($mech.stress, $mech.velocity, $mech.material, $mech.topography,
                        $ssa, $rt, $mask)
        sync!($backend)
    end
    g["membranestress! (DIVA)"] = @benchmarkable begin
        membranestress!($mech.stress, $mech.velocity, $mech.material, $mech.topography,
                        $diva, $rt, $mask)
        sync!($backend)
    end

    return g
end

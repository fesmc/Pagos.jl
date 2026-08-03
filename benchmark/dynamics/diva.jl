# ---------------------------------------------------------------------------
# The DIVA depth-integrated-viscosity chain (`src/mechanics/velocities.jl`, and
# `diva_update!` in `src/mechanics/pseudotransient.jl`).
#
# What separates DIVA from SSA in cost: the F₁/F₂ column integrals, the depth-average
# reduction behind µ̄, the β_eff correction they feed, and the 3D velocity reconstruction.
# Every one of them is a *column* operation launched one thread per column
# (`rt.launch2d` over a serial `k`-loop), not one thread per cell — a launch shape that
# behaves very differently on CPU and GPU, which is the main reason these are worth
# tracking separately from the 2D kernels.
# ---------------------------------------------------------------------------

"""
    diva_suite(fx)

Benchmarks for the viscosity integrals, the depth-average reduction, the effective-friction
correction, the full `diva_update!` chain, and the 3D velocity reconstruction.
"""
function diva_suite(fx)
    (; grid, rt, mech, mask, backend, T) = fx
    g = BenchmarkGroup()

    F1 = mech.material.viscosity_integral_1
    F2 = mech.material.viscosity_integral_2

    g["viscosity_integrals!"] = @benchmarkable begin
        viscosity_integrals!($F1, $F2, $mech, $rt, $mask)
        sync!($backend)
    end
    g["depthaverage!"] = @benchmarkable begin
        depthaverage!($mech.material.viscosity_depthaveraged, $mech.material.viscosity,
                      $rt, $mask)
        sync!($backend)
    end
    g["beta_eff_diva!"] = @benchmarkable begin
        beta_eff_diva!($mech, $rt, $mask)
        sync!($backend)
    end

    # The viscosity continuation `update_viscosity!` applies, and the `diva_update!` chain
    # that wraps it together with the integrals and β_eff.
    #
    # Both entries are configured with a `DIVAViscosityContinuation`, which is *not* the
    # solver default (`NoViscosityContinuation`, under which `diva_update!`'s first step is
    # a no-op and the chain costs only the integrals plus β_eff). The expensive branch is
    # the one worth tracking: it dominates `diva_update!` by an order of magnitude — a
    # per-layer effective strain rate, a `pow`-bound Glen kernel over the full column, and
    # a depth-average reduction — so a regression in it would be invisible under the
    # default. The solver entries in `pseudotransient.jl` run the default, so both
    # configurations are covered.
    #
    # `strainrate_reg` must be non-zero: the continuation divides by the effective strain
    # rate, which the uniform slab drives to zero.
    vc = DIVAViscosityContinuation(T; strainrate_reg = T(1.0e-6))
    solver = PseudoTransientSolver(grid; maxiter = 1, viscosity_continuation = vc)

    g["update_viscosity! (DIVA)"] = @benchmarkable begin
        update_viscosity!($mech, $vc, $rt, $mask)
        sync!($backend)
    end
    g["diva_update! (DIVA continuation)"] = @benchmarkable begin
        diva_update!($mech, $solver, $rt, $mask)
        sync!($backend)
    end

    # The Eq. 17 reconstruction: run once after a converged solve, not per iteration, so
    # its cost is amortised — but it is the only place the full 3D velocity is written, and
    # it is a serial-in-`k` column sweep.
    g["velocities3D! (DIVA)"] = @benchmarkable begin
        velocities3D!($mech, $rt, DIVAMomentumBalance(), $mask)
        sync!($backend)
    end
    g["velocities3D! (SSA)"] = @benchmarkable begin
        velocities3D!($mech, $rt, SSAMomentumBalance(), $mask)
        sync!($backend)
    end

    return g
end

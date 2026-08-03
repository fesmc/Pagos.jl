# ---------------------------------------------------------------------------
# The pseudo-transient momentum solver (`src/mechanics/pseudotransient.jl`).
#
# Two kinds of entry, deliberately kept apart:
#
#   "iteration"  — a fixed number of PT iterations (`abstol = -1`, see below). Measures
#                  **cost per iteration** and nothing else. A regression here means a
#                  kernel on the residual path got slower, full stop.
#   "solve"      — a real solve, run to `abstol`. Measures **cost to an answer**, which is
#                  per-iteration cost × iterations-to-converge. This one moves whenever the
#                  tuning, the pseudo-time-step rule or the convergence criterion changes,
#                  *without any kernel getting slower* — so read a regression here against
#                  the iteration entries before concluding anything about throughput.
#
# The fixed-iteration entries use `abstol = -one(T)`, not `0.0`: `pseudo_transient!` loops
# `while err > abstol`, and `err` is a norm, so a negative tolerance can never be met and
# exactly `maxiter` iterations run on every sample. With `abstol = 0.0` a slab that reaches
# its fixed point to the last bit — which the uniform fixture can — would exit early and
# silently turn a "fixed 25 iterations" entry into a variable one.
#
# Every solve entry runs `restore_momentum!` then `reset_velocity!` in `setup`, and
# `evals = 1` so that setup precedes every timed body rather than every *batch* of them.
# Both halves are load-bearing. `reset_velocity!` fixes the initial guess: `pseudo_transient!`
# iterates the velocity in place, so without it sample n would start from sample n-1's
# converged answer and time a solve that had already happened. `restore_momentum!` fixes
# everything else the solve reads — β_eff and the viscosities, which other benchmarks in
# this suite overwrite (see `_momentum_inputs` in `common.jl`). Without it the iteration
# count, and therefore this entry's whole meaning, would depend on the order BenchmarkTools
# happened to walk the group in.
# ---------------------------------------------------------------------------

"""
    pseudotransient_suite(fx)

Benchmarks for the PT solver: the per-iteration pieces (`pseudo_dt!`, `pseudo_rate!`),
fixed-iteration throughput for SSA and DIVA, and full solves to convergence.
"""
function pseudotransient_suite(fx)
    (; grid, rt, mech, cst, mask, backend, T) = fx
    g = BenchmarkGroup()

    ssa  = SSAMomentumBalance()
    diva = DIVAMomentumBalance()

    # ---- per-iteration pieces ------------------------------------------------
    ref_solver = PseudoTransientSolver(grid)

    # The Gershgorin bound. Not per-iteration: the loop recomputes it only when DIVA's
    # `div_update` cadence fires (a fresh β_eff would otherwise leave the bound optimistic
    # and the explicit iteration divergent), and once before the loop otherwise. Tracked on
    # its own because it is the one kernel whose cost does *not* scale with the iteration
    # count, so a regression in it hides inside a solve entry.
    g["pseudo_dt!"] = @benchmarkable begin
        pseudo_dt!($ref_solver, $mech, $cst, $rt, $mask)
        sync!($backend)
    end

    # One full residual evaluation: gradients → strain rate → membrane stress → basal
    # stress → dotvel. The bulk of an iteration's cost, minus the velocity update and the
    # convergence reduction.
    g["pseudo_rate! (SSA)"] = @benchmarkable begin
        pseudo_rate!($mech, $cst, $rt, $ssa, $ref_solver, $mask)
        sync!($backend)
    end
    g["pseudo_rate! (DIVA)"] = @benchmarkable begin
        pseudo_rate!($mech, $cst, $rt, $diva, $ref_solver, $mask)
        sync!($backend)
    end

    # ---- fixed-iteration throughput -----------------------------------------
    iters = 25
    it = g["iteration"] = BenchmarkGroup()
    fixed_ssa = PseudoTransientSolver(grid; maxiter = iters, abstol = -one(T))
    fixed_diva = PseudoTransientSolver(grid; maxiter = iters, abstol = -one(T),
                                       div_update = PeriodicDIVUpdate(10))

    it["SSA ($iters iters)"] = @benchmarkable begin
        pseudo_transient!($mech, $cst, $fixed_ssa, $rt, $ssa, $mask)
        sync!($backend)
    end setup = (restore_momentum!($fx); reset_velocity!($fx)) evals = 1 samples = 10 seconds = 20

    # DIVA under a `PeriodicDIVUpdate`: the depth-integrated-viscosity chain refreshes
    # every 10 iterations, which is the cadence a real DIVA run pays. Under the default
    # `NoDIVUpdate` the loop would never call `diva_update!` and the entry would differ
    # from the SSA one only by the per-layer effective strain rate.
    it["DIVA ($iters iters)"] = @benchmarkable begin
        pseudo_transient!($mech, $cst, $fixed_diva, $rt, $diva, $mask)
        sync!($backend)
    end setup = (restore_momentum!($fx); reset_velocity!($fx)) evals = 1 samples = 10 seconds = 20

    # ---- solves to convergence ----------------------------------------------
    # These are the most expensive entries in the suite (the DIVA one is seconds, not
    # milliseconds: `diva_update!` rewrites β_eff by ~30× mid-solve at the `div_update`
    # cadence, and the iteration has to re-converge from there). Sample counts are set
    # explicitly rather than left to the one-second default, which would take a single
    # sample of each — a minimum over one is not a minimum.
    solve_ssa = PseudoTransientSolver(grid; abstol = 1.0e-8, maxiter = 500)
    solve_diva = PseudoTransientSolver(grid; abstol = 1.0e-8, maxiter = 500,
                                       div_update = PeriodicDIVUpdate(10))
    # The autotuner is a distinct code path (a spectral estimate and a γ update folded into
    # every iteration), and the thing it is *for* is cutting the iteration count — which
    # only a converged entry can see.
    solve_auto = PseudoTransientSolver(grid; abstol = 1.0e-8, maxiter = 500,
                                       tuning = AutotunedDynamicRelaxation())

    sv = g["solve"] = BenchmarkGroup()

    sv["SSA"] = @benchmarkable begin
        pseudo_transient!($mech, $cst, $solve_ssa, $rt, $ssa, $mask)
        sync!($backend)
    end setup = (restore_momentum!($fx); reset_velocity!($fx)) evals = 1 samples = 5 seconds = 20

    sv["DIVA"] = @benchmarkable begin
        pseudo_transient!($mech, $cst, $solve_diva, $rt, $diva, $mask)
        sync!($backend)
    end setup = (restore_momentum!($fx); reset_velocity!($fx)) evals = 1 samples = 5 seconds = 20

    sv["SSA autotuned"] = @benchmarkable begin
        pseudo_transient!($mech, $cst, $solve_auto, $rt, $ssa, $mask)
        sync!($backend)
    end setup = (restore_momentum!($fx); reset_velocity!($fx)) evals = 1 samples = 5 seconds = 20

    return g
end

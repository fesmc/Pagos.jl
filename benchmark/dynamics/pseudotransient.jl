# ---------------------------------------------------------------------------
# The pseudo-transient momentum solver (`src/mechanics/pseudotransient.jl`).
#
# Two kinds of entry, kept apart:
#   "iteration" — fixed PT iterations (`abstol = -one(T)`, unreachable so `maxiter`
#                 always runs). Measures cost per iteration only.
#   "solve"     — a real solve to `abstol`. Measures cost to an answer, which moves
#                 with tuning/step-rule/convergence-criterion changes even when no
#                 kernel got slower — read a regression here against the iteration
#                 entries before concluding anything about throughput.
#
# Every solve/iteration entry runs `restore_momentum!` + `reset_velocity!` in `setup`
# (evals = 1, so setup precedes every sample): otherwise sample n would start from
# sample n-1's converged velocity/β_eff/viscosities (which other benchmarks in this
# suite overwrite, see `_momentum_inputs` in `common.jl`), making the iteration count
# depend on the order BenchmarkTools walks the group in.
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

    # `pseudo_rate!` reads the membrane prefactor cache rather than rebuilding it (see
    # `membrane_prefactors!`), and `pseudo_transient!` — not `pseudo_rate!` — is what
    # normally fills it. Building it here is not tidiness: on a fresh solver the cache is
    # zero, which makes the membrane kernel's work vanish and would time an iteration doing
    # strictly less than a real one.
    membrane_prefactors!(ref_solver, mech, rt, mask)

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

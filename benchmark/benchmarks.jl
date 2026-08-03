# ---------------------------------------------------------------------------
# Pagos benchmark suite — the PkgBenchmark entry point.
#
# Defines `SUITE`, a `BenchmarkGroup` covering the functions whose performance is tracked
# across commits. This file is *not* an exploratory script: every entry here is meant to
# stay comparable over time, so adding one is a commitment to keeping its inputs fixed.
# One-off "which primitive is faster" experiments belong in `basics/`.
#
# Run it:
#
#   julia --project=benchmark -e '
#       using PkgBenchmark, Pagos
#       r = benchmarkpkg(Pagos)
#       export_markdown("benchmark/results/HEAD.md", r)'
#
# Compare two commits (PkgBenchmark checks each one out into a temp dir itself, so the
# working tree is left alone — but commit first, it benchmarks commits, not the tree):
#
#   julia --project=benchmark -e '
#       using PkgBenchmark, Pagos
#       j = judge(Pagos, "HEAD", "main")
#       export_markdown("benchmark/results/HEAD-vs-main.md", j)'
#
# GPU entries are off unless `PAGOS_BENCH_GPU=1` and CUDA is functional — see
# `gpu_requested()` in `common.jl`.
# ---------------------------------------------------------------------------

using BenchmarkTools
using Pagos
using KernelAbstractions

# BenchmarkTools' own default is 5 s per entry, which puts this suite over eight minutes —
# too slow to run on every commit, which defeats the purpose. One second still buys ≥ 10
# samples for every kernel entry (the slowest is ~100 ms), and the few entries that need
# more say so explicitly. Set before the suite is built: `@benchmarkable` snapshots these
# defaults into each benchmark at construction.
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1.0

const SUITE = BenchmarkGroup()

include("common.jl")

include("numerics/differences.jl")
include("numerics/pointwise.jl")
include("numerics/integrators.jl")

include("dynamics/strainrate.jl")
include("dynamics/stress.jl")
include("dynamics/diva.jl")
include("dynamics/pseudotransient.jl")

"""
    build_suite!(suite, fx)

Fill `suite` with the `numerics` and `dynamics` groups built against the fixture `fx`.
Shared by the CPU suite and the gated GPU one, so the two measure the same call sites on
different backends rather than drifting apart.

!!! note "Not yet covered"
    Topography (advection, calving, masks), material laws (rate factor, creep, flow law)
    and thermodynamics are deliberately out of scope for now — the tracked suite starts at
    the momentum hot path. They are expected to join later; when they do, they want their
    own top-level group rather than being folded into `dynamics`.
"""
function build_suite!(suite, fx)
    num = suite["numerics"] = BenchmarkGroup()
    num["differences"] = differences_suite(fx)
    num["pointwise"]   = pointwise_suite(fx)
    num["integrators"] = integrators_suite(fx)

    dyn = suite["dynamics"] = BenchmarkGroup()
    dyn["strainrate"]      = strainrate_suite(fx)
    dyn["stress"]          = stress_suite(fx)
    dyn["diva"]            = diva_suite(fx)
    dyn["pseudotransient"] = pseudotransient_suite(fx)

    return suite
end

build_suite!(SUITE, bench_fixture())

if gpu_requested()
    using CUDA
    if CUDA.functional()
        SUITE["gpu"] = BenchmarkGroup()
        build_suite!(SUITE["gpu"], bench_fixture(; backend = CUDABackend()))
    else
        @warn "PAGOS_BENCH_GPU is set but CUDA.functional() == false — skipping GPU group."
    end
end

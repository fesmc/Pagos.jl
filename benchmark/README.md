# `benchmark/`

Two different jobs live here, and keeping them apart is the point.

| | question | where |
| --- | --- | --- |
| **tracked suite** | is *this* implementation getting slower? | `benchmarks.jl` + `numerics/` + `dynamics/` |
| **design experiments** | which of these alternatives should we build on? | `basics/` (see its own README) |

## The tracked suite

`benchmarks.jl` defines `SUITE`, a `BenchmarkTools.BenchmarkGroup`, in the layout
[PkgBenchmark](https://github.com/JuliaCI/PkgBenchmark.jl) expects. Everything it times is
built from closed-form inputs at one fixed size, so two runs are comparable if and only if
the code changed.

```
SUITE
├── numerics
│   ├── differences    ∂x!, ∂y!, ∂x₁₂!, ∂x₃! (uniform and sigma)
│   ├── pointwise      Pagos.pointwise!, against a broadcast reference
│   └── integrators    step! per fixed-step scheme
├── dynamics
│   ├── strainrate     velocity gradients, raw/effective strain rates, membrane stress
│   ├── stress         driving, basal and 3D deviatoric stress
│   ├── diva           F₁/F₂ integrals, depth-average, β_eff, diva_update!, velocities3D!
│   └── pseudotransient
│       ├── pseudo_dt!, pseudo_rate!
│       ├── iteration/  fixed 25 PT iterations — cost *per iteration*
│       └── solve/      run to abstol — cost *to an answer*
└── gpu                the same two groups on a CUDA backend (gated, see below)
```

### Run it

```bash
julia --project=benchmark -e '
    using PkgBenchmark, Pagos
    r = benchmarkpkg(Pagos)
    export_markdown("benchmark/results/HEAD.md", r)'
```

Roughly 2 minutes of benchmarking, plus a one-off tuning pass (~1.5 min) cached in
`benchmark/tune.json`. Both `results/` and `tune.json` are gitignored.

### Compare two commits

```bash
julia --project=benchmark -e '
    using PkgBenchmark, Pagos
    j = judge(Pagos, "HEAD", "main")
    export_markdown("benchmark/results/HEAD-vs-main.md", j)'
```

PkgBenchmark checks each commit out into a temporary directory itself, so your working tree
is untouched — but it benchmarks *commits*, so commit before comparing. To compare against
a commit from before this suite existed, hand it the current script:

```julia
judge(Pagos, "HEAD", "some-old-sha"; script = "benchmark/benchmarks.jl")
```

### GPU

Off by default, so the CPU baseline stays comparable across machines:

```bash
PAGOS_BENCH_GPU=1 julia --project=benchmark -e 'using PkgBenchmark, Pagos; benchmarkpkg(Pagos)'
```

It builds a second fixture on `CUDABackend()` and runs the *same* call sites through it.
Every timed body ends in `sync!(backend)`; without it a GPU entry would measure launch
overhead rather than execution.

## Reading the results

**`iteration/` vs `solve/` is the distinction to internalise.** `iteration/` runs a fixed
25 PT iterations, so it moves only when a kernel on the residual path gets slower.
`solve/` runs to `abstol`, so it is per-iteration cost × iterations-to-converge, and it
moves whenever the tuning, the pseudo-time-step rule or the convergence criterion changes —
with every kernel exactly as fast as before. A regression in `solve/` with `iteration/`
flat is a convergence change, not a throughput one. Look at both before drawing a
conclusion.

**Allocation counts matter as much as time** for the kernel entries. The launches should be
allocation-free bar a few hundred bytes of launch overhead; a jump means something started
boxing or capturing (`Pagos.pointwise!`'s docstring records exactly this failure mode,
which is a hard error on GPU and merely slow on CPU).

## Adding a benchmark

Adding an entry is a commitment to keeping its inputs fixed, forever — that is what makes
the history comparable. So:

- **Never `rand`.** Closed-form inputs only. Several kernels branch on their data
  (`node_active`, regularisation floors, `H > 0` guards), so random inputs let the fraction
  of cells taking each branch drift between runs.
- **Mutated state must be restored in `setup`, with `evals = 1`.** Benchmarks in a group
  run in dictionary order, which is neither source order nor stable, and the fixture is
  shared — so anything reading state another entry writes needs `restore_momentum!` (see
  `_momentum_inputs` in `common.jl`). Without it, an entry's meaning depends on what
  happened to run before it.
- **One size.** Resolution sweeps are a `basics/` question.
- **End the timed body with `sync!(backend)`** so the entry is honest on GPU.

Topography, material laws and thermodynamics are not covered yet; when they are, they want
their own top-level group rather than being folded into `dynamics`.

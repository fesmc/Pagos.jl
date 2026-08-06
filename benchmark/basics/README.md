# `basics/` — design experiments, kept as a record

These scripts answered "which primitive should this be built on?" *before* the code above
them existed. They are kept because the decisions they justify are still load-bearing —
not because they are expected to run.

| Script | Question it answered |
| --- | --- |
| `map.jl` | `map!` vs. broadcast vs. a hand-written loop vs. a KernelAbstractions kernel, including masked (`view`) variants |
| `map_ka.jl` | `map!` vs. KA for a regularized-Coulomb friction map, CPU and GPU |
| `differences.jl` | sequential vs. fused planar derivatives; cost of each `AbstractIndexing` strategy |
| `difference_kernels.jl` | the same comparison at the kernel level |
| `stress.jl` | the fused `deviatoric_stress!` kernel across 2D and 3D grids |
| `chmy_comparison.jl` | Pagos' own launcher vs. Chmy's, before the Chmy migration |
| `velocities.jl` | COO-fill + `sparse()` vs. pre-built CSC + direct `nzval` writes; UMFPACK vs. CUDSS; PT iteration vs. a direct solve |
| `pseudotransient/` | PT vs. linear solve, single- and multi-resolution, the PT loop itself |
| `gpu/` | where the DIVA/SSA PT loop's GPU time goes: launch worksize, field layout, `hlerp` cost, the per-iteration copies (see its own README) |

## They do not run as-is

**`gpu/` is the exception** — it targets the current API and is meant to be re-run whenever
the momentum path changes. Everything below applies to the scripts in this directory only.


Most of them target the pre-Chmy API (`Domain`, `State`, `IceSheet`, `RegularGrid`,
`LinearMomentumSolver2D`), which no longer exists in that form. Reviving one means porting
it, and if it is worth porting it is probably worth promoting into the tracked suite
instead.

They also use dependencies the benchmark environment no longer carries — `Chairmarks`
(`@b`), `CairoMakie`, `CUDSS`, `NCDatasets`. Add them back to `benchmark/Project.toml`, or
run the script under its own temporary environment, if you do resurrect one.

## What replaced them

`benchmark/benchmarks.jl` and the `numerics/` + `dynamics/` groups beside it: a
`BenchmarkGroup` with fixed inputs and fixed sizes, run through PkgBenchmark so
`judge(Pagos, "HEAD", "main")` reports per-function regressions across commits. Different
job — these scripts compare *alternatives at one point in time*, the suite tracks *one
implementation over time* — so both are worth having.

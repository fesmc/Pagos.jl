# `basics/gpu/` — where the DIVA/SSA PT loop's GPU time goes

Design experiments, in the sense the parent `basics/README.md` means: they compare
*alternatives at one point in time* rather than tracking one implementation over time.
Unlike the older scripts one directory up, **these run against the current API** and are
meant to be re-run when the momentum path changes.

| Script | Question it answers |
| --- | --- |
| `layout_and_launch.jl` | How much of a depth-integrated kernel's cost is the launch worksize and the field layout, rather than the kernel body? |
| `kernel_variants.jl` | With the worksize fixed, which kernel bodies on the residual path still do avoidable work? |
| `pt_loop.jl` | Stacked on the real loop, what does it add up to — and is the answer unchanged? |

```bash
julia --project=benchmark benchmark/basics/gpu/pt_loop.jl          # Float64
julia --project=benchmark benchmark/basics/gpu/pt_loop.jl f32 5    # Float32, 5 reps
```

They need a functional CUDA device. The geometry is the tracked suite's closed-form ice
slab (`benchmark/common.jl`) at 760×760×2 — no data file, and the kernels' cost does not
depend on the values, only on which `node_active` branch each node takes.

## What they found

Measured on an RTX 2070 Super Max-Q (8 GB, ~280 GB/s achievable), Julia 1.12, CUDA.jl 5.11,
Float64, 100 fixed PT iterations. `pt_loop.jl` checks every configuration against
`pseudo_transient!`'s own output: all of them are **bit-identical**.

| configuration | ms/iter | speedup |
| --- | --- | --- |
| baseline | 4.77 | 1.00× |
| + flat z-sweep | 2.40 | 1.99× |
| + no per-launch sync | 1.82 | 2.62× |
| + fused kernels & hoisted prefactors | 1.33 | **3.59×** |

Float32: 2.34 → 0.91 ms/iter (2.58×). On real AIS geometry at 8 km (761×761, roughly half
the nodes inactive) the same stack gives 2.84 → 0.97 ms/iter, i.e. 2.92× in Float64 and
2.15× in Float32.

### 1. Every depth-integrated kernel does 3× the work it needs to

Chmy's `Launcher` sweeps `size(grid, Center()) .+ 2`. On `grid2d` that is
`(nx + 2, ny + 2, 3)`: the z axis has extent 1, so the rule adds a halo ring in a dimension
that has none, and every 2D kernel reads and writes **three** k-planes instead of one.

Nothing reads those planes. The 2D grid operators (`∂x`, `∂y`, `lerp`/`hlerp` at
`aa`/`ab`/`acx`/`acy`) never index `k ± 1`, and the column kernels that read a
depth-integrated field index it explicitly — `_dotvel_staggered_bp!` writes
`lerp(H, NODE_ACX, grid2d, i, j, 1)`.

`FlatLauncher` in `common.jl` is the one-object fix: same kernels, worksize
`(nx + 2, ny + 2, 1)` and `Offset(-1, -1, 0)`. Per-kernel, 2.9–4.2×; on the whole loop, 1.99×.

The same halo convention makes a `grid2d` `Field` a single k-plane of an
`(nx + 4, ny + 4, 5)` parent — **5.1× the memory its interior needs**, 22 MB per 2D field at
this size. That, not the column, is what puts a Float64 AIS state at 8 km over an 8 GB card.

### 2. Chmy synchronizes after every launch

`Chmy/src/KernelLaunch.jl:117` calls `KernelAbstractions.synchronize(backend)` at the end of
every `Launcher` call, so no two kernels ever overlap and each launch's host cost is fully
exposed (~5 launcher calls per PT iteration, ~9 µs per `cuLaunchKernel`). Dropping it is
worth a further 1.3×. It is a Chmy policy question, not a Pagos one — and the gain is
size-dependent: at 381×381 the async queueing costs more than it saves.

### 3. `hlerp` is nine FP64 divisions a node, recomputed every iteration

`_membrane_stress_staggered!` is the most expensive kernel in the Float64 loop — 32% of
device time before any of this. It is not bandwidth: at the flat worksize it runs at ~40%
of the roofline while `_dotvel_staggered!` runs at ~160 GB/s.

Chmy's `HarmonicLinear` rule is `inv(muladd(t, inv(b) - inv(a), inv(a)))`, applied three
times for a 2D interpolation — nine divisions, at FP64's 1/32 rate on a consumer card. The
Float64/Float32 cost ratio gives it away: 3.8× for this kernel where bandwidth alone would
predict 2.0× (and 9.3× for `gershgorin_dt!`, which is worse but only fires every 50
iterations).

Collapsing `hlerp` to the 4-point harmonic mean `4 / Σ 1/ηᵢ` — five divisions instead of
nine — buys almost nothing (1.05×). Hoisting the prefactor out of the loop buys **4.0×**
in Float64 against only 1.28× in Float32, which is the same story from the other side:
both `2ηH` at `aa` and `hlerp(η)·lerp(H)` at `ab` depend only on `η` and `H`, and neither
moves during a DIVA solve (`_iterate_viscosity!` is a no-op for DIVA; `NoDIVUpdate` leaves
`µ̄` alone). The build costs 0.7 of an iteration and a 200-iteration solve recomputes it 200
times for one answer.

### 4. Four `copyto!`s a iteration, at 60% of copy bandwidth

- `update_basalstress!` copies `ū` into `velocity.base_{x,y}` and reads it straight back.
  Under the SSA limit `u_b = ū` the kernel can read `ū` itself and write `base_{x,y}` on the
  way past: 1.33×.
- `copyto!(u_old, u)` then `u = u_old + θ·dv·dτ` reads `u` twice. One kernel doing both, for
  both components at once: 1.79×.

Separately, `asarray(f) = interior(f)` is a **strided** `SubArray`, so every `copyto!` and
every broadcast on it goes through a Cartesian-index kernel: 156 GB/s against 243 GB/s for
the same bytes as the contiguous parent k-plane `view(parent(f), :, :, 3)`. That is a 1.6×
lever on any `asarray` op that survives the fusions above.

## Two measurement traps

**Warm up past the tuning cadence.** `AutotunedDynamicRelaxation(cadence = 50)` first calls
`_arm_tuning`/`_tune!` at iteration 50, so a 5-iteration warm-up leaves them and their
three-array `mapreduce`s uncompiled. Timing 100 iterations after a 5-iteration warm-up gives
44.4 ms/iter against 2.45 ms/iter fully warm — ~4 s of GPUCompiler inside the timed solve.
`docs/src/examples/ais-momentum/cpu-gpu.jl` warms up with `maxiter = 5` and hits exactly
this; its published GPU number is ~8× too slow.

**`set_theme!(theme_latexfonts())` breaks CUDA initialisation.** CairoMakie's font setup
leaves the process in a state where a later CUDA driver init fails with
`CUDA_ERROR_NOT_INITIALIZED`. `docs/src/examples/ais-momentum/helpers.jl` ends with that
call and `cpu-gpu.jl` does `using CUDA` *after* including it, so `CUDA.functional()` is
`false` and the whole GPU comparison is silently skipped. Touch the device before the
`include`, or move the `set_theme!`.

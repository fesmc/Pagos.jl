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

**Landed in `src`.** `Pagos.FlatLauncher` (`src/api/runtime.jl`) is this type, hardened for
library use — it implements the `bc =` path Chmy's `Launcher` has (and this probe's version
throws on), and `Runtime` installs it as `rt.launch2d`, so all ~40 depth-integrated launch
sites got the fix without a call-site edit. The probe copy here stays as the measurement
apparatus. `outer_width` still falls back to a plain `Launcher`.

Cross-checked on **CPU** (4 threads, 380×380×2, Float64, 100 fixed PT iterations, best of 5)
by swapping `rt.launch2d` back to a Chmy `Launcher`: 15.83 → 7.13 ms/iter, **2.22×**, with
`depthaverage_x`/`_y` bit-identical between the two. So this is not a GPU-only effect — the
redundant planes are redundant *work*, and both backends were doing 3× of it.

The same halo convention makes a `grid2d` `Field` a single k-plane of an
`(nx + 4, ny + 4, 5)` parent — **5.1× the memory its interior needs**, 22 MB per 2D field at
this size. That, not the column, is what puts a Float64 AIS state at 8 km over an 8 GB card.
Reclaiming it is now a Pagos-side change (per-axis `halo` at the allocation sites in
`src/api/state.jl`) rather than an upstream blocker — see `pagos-roadmap/memreduce.md` §2.

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

**Landed in `src`.** The prefactor pair lives on the solver
(`PseudoTransientSolver.membrane_pre_aa`/`membrane_pre_ab`), is written by
`membrane_prefactors!` and read by a cached `membranestress!` method
(`src/mechanics/strainrate.jl`). Rather than assuming `µ̄` is static, the cache is
*invalidated where `µ̄` is written*: `_refresh_membrane_prefactors!` runs after
`_iterate_viscosity!` and after `diva_update!`, and dispatches to a no-op under
`NoViscosityContinuation`. That makes it correct for SSA-with-continuation too, not just the
DIVA/`NoDIVUpdate` case this section measured.

Whole-loop on **CPU** (4 threads, 380×380×4, Float64, 100 fixed iterations, best of 5),
verified bit-identical to the old path in all four configurations:

| configuration | before | after | |
|---|---:|---:|---:|
| DIVA + `NoDIVUpdate` | 10.93 | 7.84 ms/iter | **1.39×** |
| SSA + `NoViscosityContinuation` | 9.52 | 7.79 ms/iter | **1.22×** |
| DIVA + `PeriodicDIVUpdate(10)` | 15.29 | 13.11 ms/iter | 1.17× |
| SSA + `GlenViscosityContinuation` | 15.01 | 13.81 ms/iter | 1.09× |

Smaller than the 4.0× above because that figure is this *kernel* in isolation on a GPU,
where FP64 division runs at 1/32 rate; on CPU the divisions are cheaper and the kernel is one
of ~7 in the loop. The bottom two rows rebuild the cache as often as `µ̄` moves, so they gain
only from the work splitting into two simpler kernels — but they do gain, so the cache costs
nothing anywhere. The GPU figure has not been re-measured against `src` (no working CUDA on
the dev box at the time).

### 4. Four `copyto!`s a iteration, at 60% of copy bandwidth

- `update_basalstress!` copies `ū` into `velocity.base_{x,y}` and reads it straight back.
  Under the SSA limit `u_b = ū` the kernel can read `ū` itself and write `base_{x,y}` on the
  way past: 1.33×. **In `src`**: `update_basalstress!`'s `ActiveFrictionUpdate` method
  (`src/mechanics/pseudotransient.jl`) is one `_basalstress_active!` launch now, reading
  `velocity.depthaverage_{x,y}` once and writing both `velocity.base_{x,y}` and
  `stress.base_{x,y}` from it — no `copyto!`, no separate `basalstress!` call.
- `copyto!(u_old, u)` then `u = u_old + θ·dv·dτ` reads `u` twice. One kernel doing both, for
  both components at once: 1.79×. **In `src`**: `pseudo_transient!`'s `MomentumBalance2D`
  loop calls `_pseudo_vel_and_store!` (`src/mechanics/pseudotransient.jl`), which reads each
  component's pre-update value once and writes both `u_old` and the relaxed `u` from it. The
  `MomentumBalance3D` (Blatter-Pattyn) loop still uses the unfused `copyto!` + `pseudo_vel!`
  pair — this fusion was only ever validated (bit-for-bit, against `pseudo_transient!`
  itself, on real 760×760 GPU geometry — see `pt_loop.jl`) for the 2D SSA/DIVA path.

  !!! warning "The isolated 1.33×/1.79× numbers above assume the flat worksize"
      `kernel_variants.jl` measures both fusions at `FlatLauncher`'s `(nx+2, ny+2, 1)`
      worksize, not through `rt.launch2d`. `pseudo_vel!`'s `copyto!`/broadcast pair paid no
      `Launcher` tax at all (`asarray`/`interior` sweeps exactly the field's real extent) —
      routing the fused replacement through the ordinary `rt.launch2d` instead (§1's
      `(nx+2, ny+2, 3)`) makes it pay a 3× tax on work that previously paid none, which
      **measured ~10% slower** for the whole PT loop, not faster. Both `src` methods above
      launched by hand at the flat worksize instead. **Superseded**: finding #1's fix has
      since landed — `rt.launch2d` is a `Pagos.FlatLauncher` (`src/api/runtime.jl`) with
      exactly that worksize, so both methods now launch through it like every other kernel
      and the hand-launch helper is gone.
      With that fix, the real `pseudo_transient!` loop measures 4.09 → 3.76 ms/iter (median
      of 15, Float64, same 760×760 slab and RTX 2070 Super Max-Q as above) from these two
      fusions alone — smaller than the isolated numbers since they are 2 of ~7 kernels in
      the unflattened loop, but a real, reproducible ~1.09× on top of the unchanged
      3-k-plane launcher, sync-per-launch and `hlerp`-recompute costs findings #1–#3
      describe.

Separately, `asarray(f) = interior(f)` is a **strided** `SubArray`, so every `copyto!` and
every broadcast on it goes through a Cartesian-index kernel: 156 GB/s against 243 GB/s for
the same bytes as the contiguous parent k-plane `view(parent(f), :, :, 3)`. That is a 1.6×
lever on any `asarray` op that survives the fusions above — moot for the two above once
fused into kernels (a `Chmy.Field`'s own `getindex`/`setindex!` isn't the strided path this
measures), so there was nothing left in the per-iteration loop to apply it to. The one
`asarray`-based `copyto!` this doesn't reach is `velocities3D!`'s `surface_{x,y} ←
depthaverage_{x,y}` (`src/mechanics/velocities.jl`) — called once per solve, not once per
iteration, so out of scope for what this benchmark measured.

### 5. Groupsize: type-parameter it, but don't expect it to matter much

`FlatLauncher` computed `heuristic_groupsize(backend, Val(3))` once at construction — right —
but stored the result in a struct *field* and read it back at every launch. A field load
isn't a compile-time constant, so `StaticSize(launcher.groupsize)` couldn't fold, and the
launch config never reached the type system on `rt.launch2d` — the ~40 depth-integrated
launch sites, i.e. every per-iteration kernel in the PT loop. `Worksize` was already a type
parameter; `GroupSize` now is too (`src/api/runtime.jl`), the same specialization a plain
`Chmy.Launcher` gets for free from `heuristic_groupsize` folding on the backend's type.
**Landed in `src`.**

Given that, the natural next question — is the default block shape actually a good one for
the kernel that matters most? — turned out to have a non-answer here. Swept
`(32,8,1)`/`(16,16,1)`/`(64,4,1)`/`(128,2,1)`/`(256,1,1)`/... against
`_membrane_stress_staggered!` (§3, the single most expensive kernel) at the fixed flat
worksize; with drift controlled for (see the third measurement trap below), every candidate
lands within **1.13×** of every other — one weak outlier at `(32,32,1)` aside, groupsize is
simply not a lever on this kernel on this card. `kernel_variants.jl`'s §0 keeps the sweep
(round-robin, so a re-run on different hardware is trustworthy without redesigning it), but
the honest conclusion is "shipped `GroupSize` as a type parameter because it's free and
correct, not because a specific block shape was worth hardcoding."

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

**`bench`'s consecutive batches measure clock drift when comparing many candidates in a
row, not just two.** Fine for A/B (§1–§4: two variants, back-to-back, similar recent
thermal history). Broke down sweeping 12 groupsizes (§5): the *same* `(32, 8, 1)` config,
run at different positions in the sequence, measured anywhere from 315 us to 595 us — a
1.9× spread with the config held fixed — and this survived a confirmed thermal soak
(`nvidia-smi`: 930 → 1395 MHz, stable throughout). This card drifts under sustained load
regardless of temperature having settled; a naive first pass at the §5 sweep reported a
bogus "1.88× best-vs-heuristic" that was pure position artifact. Fix: round-robin — one
short timing per candidate per round, cycling through all candidates for many rounds — so
every candidate samples the same drift trajectory and it cancels out of the comparison
instead of aliasing onto whichever one happened to run first (or, for a two-way comparison,
interleave the two rather than timing one fully then the other).

#=

# ISMIP-HOM: shared setup

The reference-data loader, the Pagos Blatter-Pattyn harness (CPU or GPU), and the plotting
vocabulary. `exp-a.jl` `include`s this file and then does only what is specific to
experiment A — the only experiment run here; see "The experiments" below for why B is not.

The benchmark is Pattyn et al. (2008), *The Cryosphere* **2**, 95–108; the data is
supplement 2 (`tc-2007-0019-sp2`), the format specification supplement 3.

## The experiments

Both are **diagnostic** (no time evolution) and share a geometry: a parallel-sided slab of
mean thickness `H = 1000 m` on a bed of mean slope `α = 0.5°`, with a surface
`z_s(x, y) = -x \tan\alpha` and a sinusoidal bed of 500 m amplitude and wavelength `L`,
`ω = 2π/L`:

```math
z_b^{A}(x, y) = z_s - 1000 + 500\,\sin(\omega x)\,\sin(\omega y),
\qquad
z_b^{B}(x)    = z_s - 1000 + 500\,\sin(\omega x).
```

So **B is A with the `y` dependence removed** — bumps become ripples, and the experiment
becomes tractable for flowline models. That is also why **only A is run here**: B exists so
that flowline models could take part in the intercomparison, and it is a strict
simplification of A rather than an independent test. A first-order solver that handles A
handles B by construction. Both are **no-slip** (`v(z_b) = 0`, i.e. `β² → ∞`)
with periodic lateral boundary conditions, isothermal Glen rheology `n = 3`,
`A = 10⁻¹⁶ Pa⁻ⁿ a⁻¹`. Shrinking `L` at fixed `H` is what drives the higher-order response:
at `L = 160 km` the bed undulation is shallow-ice-like, at `L = 10 km` the bed slope reaches
~17° and longitudinal stress gradients dominate.

`L = 5 km` is **deliberately excluded**. Pagos is not grid-converged there at the
resolutions these scripts can afford — the `nx = 40` and `nx = 80` solves disagree by 14 %
in mean `vx`, against ≤ 0.4 % at every other length — so a number reported at `L = 5 km`
would be measuring the grid, not the physics. The 160 → 10 km sweep is grid-converged and is
what the conclusions rest on.

## The reference: a full-Stokes ensemble, plotted as a band

`ismip_all/` holds 25 submissions of varying kind. [`MODELS`](@ref) below selects the
**full-Stokes** ones — the models that solve the momentum balance without approximation.
They are the only fair yardstick for a first-order model, and the honest reading of "truth"
is not one of their curves but the **spread between them**: where they separate, no
approximate model can be scored more finely than that gap.

So the reference is drawn as a shaded min–max envelope across the ensemble rather than as
one line per model. That is also what makes the ensemble extensible — a sixth submission
widens the band instead of needing a sixth linestyle, which four-value `linestyle` vocabulary
could not have absorbed.

Building an envelope does require putting the submissions on a common abscissa, which the
per-model profile plots never did: they submit different grids (41×41 cell centres for
`rhi3`, 61×61 nodes for `oga1`). [`interp_periodic`](@ref) does that with periodic linear
interpolation, and only for the band — every printed number is computed on each model's own
abscissa.

## File format, and three ways submissions differ

Files are `NNNMELLL.txt` — three-letter author code, model number, experiment letter,
zero-padded `L` in km. Columns are whitespace-separated, velocities in m a⁻¹, stresses and
pressures in kPa, and `x̂`/`ŷ` are normalized by `L`:

| experiment | columns |
|---|---|
| **A** (8) | `x̂  ŷ  vx(zs)  vy(zs)  vz(zs)  τxz(zb)  τyz(zb)  Δp` |

(Experiment B's five-column layout, `x̂  vx  vz  τxz  Δp`, is documented in the specification
but not read here.)

`Δp = p_I − p_H` is the departure of the isotropic pressure from hydrostatic at the bed; it
is identically zero under the shallow-ice approximation, which makes it the sharpest
higher-order diagnostic in the set.

Beyond the columns, **three things vary between submissions, and all three are detected
rather than assumed** — each would silently produce a plausible-looking wrong figure:

 1. **The file prefix is not the directory name.** 13 of the 25 directories drop the model
    number that the files carry (`oga/oga1a160.txt`), and one changes case outright
    (`yuv/Yko1a160.txt`). So files are found by globbing the directory, not by building the
    name from it.
 2. **Which coordinate runs fastest.** `rhi3` writes `x` fastest, `oga1` writes `y` fastest.
    Reshaping under the wrong convention transposes the field — invisible in the bed
    geometry, which is symmetric under `x ↔ y`, but wrong in `vx`, which is not.
 3. **Grid layout and size.** `rhi3` submits 41×41 cell centres (`x̂ ∈ (0, 1)`), `oga1`
    submits 61×61 nodes including both periodic endpoints (`x̂ ∈ [0, 1]`).

!!! warning "Submissions report `Δp` with opposite signs"
    Verified in the raw files, upstream of anything these scripts do. The cleanest evidence
    is in the experiment-B files (not read by this script, but on disk) at
    `L = 5 km`: `rhi3` reports `Δp = +62.7 kPa` at `x̂ = 0` and `−63.6` at `x̂ = 0.5`; `oga1`
    reports `−63.3` and `+64.1` at the same places. The magnitudes agree to ~1 %, and
    `vx`/`vz`/`τxz` agree to ~0.3 %, so this is a reporting convention, not a disagreement
    about the physics — one of the two has `p_H − p_I` where the specification asks for
    `Δp = p_I − p_H`.

    A min–max band across models with mixed conventions would be meaningless, so
    [`dp_orientation`](@ref) detects each submission's sign (by projecting `Δp` onto the
    bed's fundamental mode, which is grid-independent) and the `Δp` band is built from
    sign-normalised references. This is **stated, not silent**: the panel says so, the
    detected signs are printed, and **Pagos is never flipped** — whether BP recovers the
    sign unaided is a result, not a plotting decision.

    On the physics, `rhi3` looks like the spec-conforming one: it puts the compressive
    (positive) anomaly on the stoss side, where the ice is being forced to climb the bed
    fastest. That is an inference from the flow direction, not something the supplements
    state, so it is what the normalisation targets and nothing more is claimed from it.
=#
using Pkg
Pkg.activate(joinpath(@__DIR__, "../../.."))   # docs/, three up from examples/ismip-hom/
using Pagos, CairoMakie, Printf, Statistics

#=
## Where the data lives

Override any of these before `include`ing an experiment script.
=#
ISMIP_DIR = @isdefined(ISMIP_DIR) ? ISMIP_DIR :
    "/home/jan/pCloudSync/PhD/Projects/Ice-Sheet-Modelling/ice-data-pagos/ismip-hom/" *
    "tc-2007-0019-sp2/ismip_all"

## Full-Stokes submissions only — see "The reference" above. The first entry is the one the
## per-model maps use as the primary. `cma1` is deliberately absent: it submits a 21×21 grid,
## the coarsest in the set and a quarter the resolution of the others, and it is the outlier
## on every diagnostic here — the widest `vx` departure at `L = 160 km` (39.42 against
## 40.16–40.48) and the worst `Δp`. A min–max band is set by its extremes, so an
## under-resolved member does not add information about the full-Stokes answer, it just
## widens the target.
MODELS = @isdefined(MODELS) ? MODELS : ["rhi3", "oga", "aas2", "rhi1"]

const LENGTHS_KM = [160, 80, 40, 20, 10]
const SLICE_Y = 0.25

figdir = joinpath(@__DIR__, "figs")
mkpath(figdir)

#=
Globbed rather than constructed: see difference (1) above. The regex is anchored and the hit
count is asserted, so a directory that breaks the naming pattern in some *fourth* way errors
here instead of silently loading a neighbouring domain length.
=#
function expfile(dir, ex, L)
    d = joinpath(ISMIP_DIR, dir)
    isdir(d) || error("no such ISMIP-HOM model directory: $d")
    pat = Regex("^[A-Za-z]+[0-9]*" * ex * @sprintf("%03d", L) * "\\.txt\$")
    hits = filter(f -> occursin(pat, f), readdir(d))
    length(hits) == 1 ||
        error("expected exactly one experiment-$ex L=$L file in $d, found $hits")
    return joinpath(d, only(hits))
end

## The submission's own identifier — the file prefix, not the directory name.
model_label(dir) = (m = match(Regex("^([A-Za-z]+[0-9]*)a160\\.txt\$"),
                              basename(expfile(dir, "a", 160))); m.captures[1])

#=
Hand-parsed rather than `DelimitedFiles.readdlm`: the files are space-aligned with variable
run lengths (and `oga1` writes Fortran `0.23760864E+02` exponent form, which `parse(Float64,
_)` handles and a naive column-width split would not), they are small enough that the
allocation pattern is irrelevant, and `DelimitedFiles` is not a dependency of the docs
environment. The column-count check is what catches a truncated or mis-named file — a
silently short row would otherwise reshape into a plausible-looking wrong picture.
=#
function read_ismip_table(path)
    isfile(path) || error("ISMIP-HOM file not found: $path")
    rows = Vector{Float64}[]
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, '#')) && continue
        push!(rows, parse.(Float64, split(s)))
    end
    isempty(rows) && error("no data rows in $path")
    ncol = length(first(rows))
    all(r -> length(r) == ncol, rows) ||
        error("ragged table in $path (expected $ncol columns on every row)")
    return reduce(vcat, permutedims.(rows))
end

#=
Experiment A is returned both as full 2D fields (for the maps) and as the `ŷ ≈ 0.25` slice
(for the profiles). The `x_fastest` test is difference (2) above: if the first coordinate
does not change between the first two rows, it is the *slow* one and the reshape has to be
transposed to reach a `[ix, iy]` field.

The slice convention follows the intercomparison's own results supplement: `sin(ω y) = 1`
exactly at `y = L/4`, so the `ŷ = 0.25` slice of experiment A cuts the bed bumps at their
**full 500 m amplitude**, which is what makes it directly comparable to experiment B's
ripples. `oga1`'s node grid contains `0.25` exactly; `rhi3`'s cell centres straddle it, so
the nearest row is used, as the intercomparison's own figures do.
=#
function load_expA(dir, L)
    path = expfile(dir, "a", L)
    tbl = read_ismip_table(path)
    nrow, ncol = size(tbl)
    ncol == 8 || error("experiment A expects 8 columns, got $ncol in $path")
    n = isqrt(nrow)
    n^2 == nrow || error("experiment A expects a square grid, got $nrow rows in $path")

    x_fastest = tbl[2, 1] != tbl[1, 1]
    fld(c) = x_fastest ? reshape(tbl[:, c], n, n) :
                         permutedims(reshape(tbl[:, c], n, n))
    x = fld(1)[:, 1]
    y = fld(2)[1, :]
    j = argmin(abs.(y .- SLICE_Y))
    F = (; vx = fld(3), vy = fld(4), vz = fld(5), txz = fld(6), tyz = fld(7), dp = fld(8))
    slice = (; x, vx = F.vx[:, j], vz = F.vz[:, j], txz = F.txz[:, j], dp = F.dp[:, j])
    return (; n, x, y, F..., slice, y_slice = y[j], x_fastest)
end

#=
## Running Pagos' Blatter-Pattyn solver on the same experiments

`roadmaps/blatter-pattyn.md` Phase 3.

!!! warning "The PT loop below is written out here rather than calling `pseudo_transient!`"
    ISMIP-HOM A and B are **periodic** domains, and Pagos cannot express that:
    `Chmy.BoundaryConditions.batch_impl` has no method for a `Periodic` axis, so `bc!` would
    `MethodError` on *every* field, and [`StaggeredGrid`](@ref) therefore rejects the topology
    at construction (see its docstring). `pseudo_transient!` in turn hardcodes
    `bc!(…, Neumann())` on the velocity every iteration, and a zero-gradient side wall is
    simply the wrong condition here — the ISMIP-HOM solution has `∂u/∂x ≠ 0` at `x = 0`.

    So the loop is written out below **using only exported building blocks**
    ([`velocitygradients!`](@ref), [`terrain_metric_correction!`](@ref),
    [`update_viscosity!`](@ref), [`membranestress!`](@ref), [`update_basalstress!`](@ref),
    [`dotvel!`](@ref), [`pseudo_vel!`](@ref), [`pseudo_dt!`](@ref)) with a periodic halo
    refresh substituted for the `bc!` call. The *discretization under test is untouched* —
    `dotvel!` assembles the same BP residual, on the same grid, with the same Gershgorin
    `Δτ`, in the same order as [`pseudo_rate!`](@ref). What is given up is
    [`AutotunedDynamicRelaxation`](@ref), whose arming/tuning hooks are private; damping here
    is a hand-set `gamma`.

    This is a workaround, not a design: the real fix is periodic halo filling implemented
    Pagos-side, which is what the `StaggeredGrid` docstring already says periodic domains need.

    **The loop is validated against the library**, not merely assumed equivalent: on a uniform
    slab — the one geometry where `Neumann` is *exact*, since the solution is `x`-independent —
    it reproduces `pseudo_transient!`'s converged profile layer-by-layer to all printed digits,
    and both match the closed form of `test/mechanics/blatterpattyn_staggered.jl`.
=#
const T_PAGOS = Float64
const SLOPE = tan(0.5 * pi / 180)     # ISMIP-HOM α = 0.5°
const A_GLEN = 1e-16                  # Pa⁻ⁿ a⁻¹, uniform (spec §2.1)
const MU_GUESS = 2e6                  # Pa a — warm-start viscosity, ≈ Glen at ε̇ ~ 5e-2 /a
const BETA_NOSLIP = 1e6               # Pa a m⁻¹; see the note below on "no slip"
const PAGOS_NZ = 12
const PAGOS_NX = Dict(:A => 32, :B => 48)

#=
### CPU or GPU — and why these runs default to threaded CPU

The GPU path works and is exact: it reproduces the CPU solve to the last printed digit
(`mean vx = 27.31`, `max = 62.90`, `err = 1.38e-01` from both, at `L = 40 km` after 2200
iterations). The state is **built on the CPU and moved in one bulk `Adapt.adapt`**, never
filled element-by-element on the device: CUDA.jl disallows scalar indexing on a `CuArray`
precisely because each write would become its own kernel launch (`ais-momentum/cpu-gpu.jl`
develops this at length). `Field` carries its own `Adapt` rule, so the whole
[`MechanicState`](@ref) moves in one call.

It is also, at ISMIP-HOM sizes, **the slow option**. Warm cost per pseudo-transient
iteration, experiment A at `L = 40 km`, measured on an RTX 2070 Super against 12 CPU threads:

| | `nx = 32` (~12k cells) | `nx = 64` (~49k cells) |
|---|---|---|
| CPU, 1 thread | 7.05 ms | 28.19 ms |
| **GPU** | 7.32 ms | 30.07 ms |
| **CPU, 12 threads** | **2.21 ms** | **11.70 ms** |

The GPU does not beat a *single* CPU thread, and threading beats it by 3.2×. These grids are
simply too small: experiment A at `nx = 32` is ~12k cells spread over ~15 kernel launches per
iteration, which is far below what it takes to fill this device. So the default here is
threaded CPU — **run these scripts with `julia -t auto`**, which is where the speed actually
comes from. Set `USE_GPU = true` before `include`ing to use the device path anyway; it is
kept working and validated because the AIS-scale problems these kernels also serve are three
orders of magnitude larger, where the ranking reverses.
=#
USE_GPU = @isdefined(USE_GPU) ? USE_GPU : false
USE_GPU && (@eval using CUDA)

Threads.nthreads() == 1 && !USE_GPU && @warn "running single-threaded — start Julia with \
    `-t auto` for a ~3x faster solve (see the CPU/GPU note in this file)"

#=
`fill_analytic!`-style filling of **interior and halo**. This is not a convenience: `setdata!`
and broadcast `.=` write the interior only, leaving the ghost ring at its allocation zero, and
`hlerp` of a zero viscosity is `NaN`, not zero — the trap [`membranestress!`](@ref) documents.
Every static input is periodic (or, for the surface, linear with a periodic *gradient*, which
is all the solver reads), so evaluating its own formula at the ghost coordinates *is* the
periodic continuation. Runs on the CPU state, before the adapt.
=#
function fill_xy!(f, g, fun)
    loc = location(f)
    ni, nj, nk = size(interior(f))
    for k in 1:nk, j in -1:(nj + 2), i in -1:(ni + 2)
        x, y, _ = coord(g, loc, i, j, k)
        f[i, j, k] = fun(x, y)
    end
    return f
end

#=
Periodic in `x`/`y`, Neumann in `z` — the same `z` condition `bc!` applies, kept so the only
difference from the library path is the horizontal one.

Written as **whole-array slice assignments on `parent(asarray(f))`** (halo width 2, so field
index `i` is parent index `i + 2`). That is what makes it backend-agnostic: the same three
lines are a handful of `copyto!`s on the CPU and a handful of kernel launches on the GPU,
with no scalar indexing on either. An `mod`-based element loop, which is what this replaced,
would have been a hard error on a `CuArray`.

The `x` pass runs over every `j` including the `y` ghosts, whose values are still garbage at
that point; the `y` pass then overwrites every `j` ghost by copying from `j` rows that the
`x` pass has already made correct, so the **corners** come out right without a third pass.

The master range is `1:nx`, `1:ny` in *field* indices. On a `Vertex` axis that leaves the
duplicate node `nx + 1` outside the master range, so it is refreshed from node `1` — correct,
because under periodicity they are the same physical degree of freedom, and it is why the
right-hand source range is sized from `ni` rather than assumed to be two columns wide.
=#
function periodic_halo!(f, nx, ny)
    a = asarray(f); p = parent(a); ni, nj, nk = size(a)
    ## x: ghosts at parent 1:2 come from field nx-1:nx, ghosts at parent nx+3:end from field 1:…
    @views p[1:2, :, :] .= p[(nx + 1):(nx + 2), :, :]
    @views p[(nx + 3):(ni + 4), :, :] .= p[3:(ni + 2 - nx + 2), :, :]
    ## y: same, now over every i, so the corners inherit the corrected x ghosts.
    @views p[:, 1:2, :] .= p[:, (ny + 1):(ny + 2), :]
    @views p[:, (ny + 3):(nj + 4), :] .= p[:, 3:(nj + 2 - ny + 2), :]
    ## z: Neumann, matching `bc!`.
    @views p[:, :, 1] .= p[:, :, 3]
    @views p[:, :, 2] .= p[:, :, 3]
    @views p[:, :, nk + 3] .= p[:, :, nk + 2]
    @views p[:, :, nk + 4] .= p[:, :, nk + 2]
    return f
end

#=
One damped pseudo-transient sweep, i.e. [`pseudo_rate!`](@ref)'s body followed by the velocity
update. Three details are load-bearing.

[`terrain_metric_correction!`](@ref) follows [`velocitygradients!`](@ref) exactly as it does
in `pseudo_rate!`. It is not optional here: it converts the horizontal gradients from
constant-`ζ` to constant-`z` derivatives, and ISMIP-HOM is precisely the regime where that
matters — at `L = 10 km` the bed slope `∂H/∂x` reaches 0.31, so the correction term is
*larger* than the term it corrects.

The viscosity halo is refreshed *between* [`update_viscosity!`](@ref) and
[`membranestress!`](@ref), which is why the body is spelled out rather than delegated.

And `copyto!(v_old, v)` must precede [`pseudo_vel!`](@ref), which computes
`v = v_old + θ·dv·Δτ` reading `v_old` rather than snapshotting it — omit it and `v_old` stays
at its initial zero, making the update `v = θ·Δτ·dv` outright, which freezes the solve the
moment the damped accumulator `dv` settles.
=#
function bp_iterate!(mech, solver, rt, momentum, cst, nx, ny; iters, gamma, theta_v,
                     glen::Bool, dt_refresh, ncheck, rtol, scale, a1 = true,
                     trace = nothing, stagnation = true)
    vx, vy = mech.velocity.x, mech.velocity.y
    err, errprev, used = Inf, Inf, 0
    for iter in 1:iters
        used = iter
        velocitygradients!(mech.velocity, mech.topography.thickness, rt)
        ## `a1 = false` is a diagnostic switch, not an option: it drops the terrain-following
        ## metric correction, which the equations doc establishes is *missing* physics. It
        ## exists so a convergence failure can be attributed to the correction or exonerated.
        if a1
            ## The correction is the first thing in this loop that reads a *gradient* at a
            ## neighbouring node — `x_dz[i, j-1]` at the `ab` nodes, `y_dz[i-1, j]` — so it is
            ## also the first thing that needs the gradients periodic across the seam. Without
            ## this the residual sticks at ~1.5 (of the driving-stress scale) on the duplicate
            ## vertex planes `i = nx+1`, `j = ny+1` while the interior converges normally,
            ## which looks exactly like a diverging solve and is not one.
            periodic_halo!(mech.velocity.x_dz, nx, ny)
            periodic_halo!(mech.velocity.y_dz, nx, ny)
            terrain_metric_correction!(mech.velocity, mech.topography.thickness,
                                       mech.topography.surface, rt)
            ## `membranestress!` then differences the corrected gradients, so those need the
            ## seam too — the correction only writes its own launch range.
            periodic_halo!(mech.velocity.x_dx, nx, ny)
            periodic_halo!(mech.velocity.y_dy, nx, ny)
            periodic_halo!(mech.velocity.x_dy, nx, ny)
            periodic_halo!(mech.velocity.y_dx, nx, ny)
        end
        if glen
            update_viscosity!(mech, solver.viscosity_continuation, rt)
            periodic_halo!(mech.material.viscosity, nx, ny)
        end
        membranestress!(mech, momentum, rt)
        update_basalstress!(mech, solver.friction_update, rt, momentum)
        dotvel!(solver.velocity_x_dt, solver.velocity_y_dt,
                mech.stress.xx, mech.stress.xy, mech.stress.xz,
                mech.stress.yy, mech.stress.yz,
                mech.stress.base_x, mech.stress.base_y,
                mech.stress.driving_x, mech.stress.driving_y,
                mech.topography.thickness, cst.density_ice, rt, momentum;
                gamma, resid_x = solver.residual_x, resid_y = solver.residual_y)
        copyto!(asarray(solver.velocity_x_old), asarray(vx))
        copyto!(asarray(solver.velocity_y_old), asarray(vy))
        pseudo_vel!(asarray(vx), asarray(solver.velocity_x_old),
                    asarray(solver.velocity_x_dt), asarray(solver.dtau_x), theta_v)
        pseudo_vel!(asarray(vy), asarray(solver.velocity_y_old),
                    asarray(solver.velocity_y_dt), asarray(solver.dtau_y), theta_v)
        periodic_halo!(vx, nx, ny); periodic_halo!(vy, nx, ny)
        glen && iter % dt_refresh == 0 && pseudo_dt!(solver, mech, cst, rt, momentum)
        if iter % ncheck == 0
            err = max(maximum(abs, interior(solver.residual_x)),
                      maximum(abs, interior(solver.residual_y))) / scale
            trace === nothing || push!(trace, (iter, err))
            (isfinite(err) && err > rtol) || break
            ## Stagnation, not convergence: the frozen-bed `β/Δζ₁` term makes `Δτ` at `k = 1`
            ## small enough that the bottom row stops moving perceptibly. Reported, not hidden.
            (stagnation && abs(err - errprev) < 1e-4 * err && iter > 4000) && break
            errprev = err
        end
    end
    return err, used
end

#=
!!! note "`β = 1e6 Pa a m⁻¹` stands in for the spec's `β² → ∞`"
    ISMIP-HOM A/B are no-slip. Pagos applies the bed condition as a *flux*,
    `µ ∂u/∂z|_b = β u_b`, so no-slip is the `β → ∞` limit rather than a Dirichlet row, and a
    finite `β` is what the discretization can express. At `1e6` the residual slip is
    `τ_b/β ≈ 0.08 m/a` against surface speeds of 11–108 m/a, i.e. below 1 % everywhere and
    smallest exactly where the reference velocities are largest. Raising it is not free:
    `β` enters the explicit Gershgorin bound as `β/(Δζ₁ H)`, which is already the stiffest
    term in the spectrum, so `1e7` buys a factor-10-smaller slip for a factor-10-smaller `Δτ`.
=#
function run_pagos_bp(experiment, L_km; nx = PAGOS_NX[experiment], nz = PAGOS_NZ,
                      beta0 = BETA_NOSLIP, warmup = 2000, iters = 60000, gamma = 0.6,
                      theta_v = 0.6, theta_mu = 0.05, reg = 1e-8, dt_refresh = 20,
                      ncheck = 250, rtol = 1e-5, a1 = true, cfl = 0.9,
                      trace = nothing, stagnation = true)
    L = L_km * 1e3; ω = 2π / L
    ny = experiment === :A ? nx : 4          # exp B has no y dependence
    dx = L / nx
    layering = CorrectedVerticalLayering(T_PAGOS, QuadraticSigmaTransform(T_PAGOS, nz))

    ## Build and fill on the CPU — see "CPU or GPU" above — then move the whole state at once.
    grid_cpu = StaggeredGrid(T_PAGOS, L, ny * dx, dx, dx, layering)
    mech_cpu = MechanicState(grid_cpu)
    bump = experiment === :A ? (x, y) -> sin(ω * x) * sin(ω * y) : (x, y) -> sin(ω * x)
    fill_xy!(mech_cpu.topography.thickness, grid_cpu.grid2d, (x, y) -> 1000 - 500 * bump(x, y))
    fill_xy!(mech_cpu.topography.surface,   grid_cpu.grid2d, (x, y) -> -x * SLOPE)
    fill_xy!(mech_cpu.friction.beta_eff,    grid_cpu.grid2d, (x, y) -> beta0)
    fill_xy!(mech_cpu.material.rate_factor, grid_cpu.grid,   (x, y) -> A_GLEN)
    fill_xy!(mech_cpu.material.viscosity,   grid_cpu.grid,   (x, y) -> MU_GUESS)
    setdata!(mech_cpu.velocity.x, zero(T_PAGOS))
    setdata!(mech_cpu.velocity.y, zero(T_PAGOS))

    grid = USE_GPU ?
        StaggeredGrid(Arch(CUDABackend()), T_PAGOS, L, ny * dx, dx, dx, layering) : grid_cpu
    mech = USE_GPU ? Pagos.Adapt.adapt(CuArray, mech_cpu) : mech_cpu
    rt = Runtime(grid); cst = Constants{T_PAGOS}()
    momentum = BlatterPattynMomentumBalance()
    periodic_halo!(mech.velocity.x, nx, ny); periodic_halo!(mech.velocity.y, nx, ny)

    solver = PseudoTransientSolver(grid, momentum;
        viscosity_continuation = BPViscosityContinuation(T_PAGOS; n_glen = 3, theta_mu,
                                                         strainrate_reg = reg),
        friction_update = ActiveFrictionUpdate(),
        pseudo_timestep = GershgorinPseudoTimeStep(; cfl))

    drivingstress!(mech, cst, rt, momentum)
    ## Normalizer for the reported residual: the rate the driving stress alone would produce,
    ## i.e. `ScaledResidual`'s scale. The first iterate (`u = 0`) has exactly this residual.
    scale = maximum(abs, interior(mech.stress.driving_x)) / cst.density_ice
    pseudo_dt!(solver, mech, cst, rt, momentum)

    t0 = time()
    ## Two stages, because a cold start is a trap: at `u = 0` the effective strain rate is the
    ## regularization floor, so Glen returns `µ ~ 2e9 Pa a`, and a nearly rigid slab barely
    ## deforms — a self-consistent stiff state the continuation crawls out of. Stage 1 holds
    ## `µ` at `MU_GUESS` to get a velocity field of the right order first.
    bp_iterate!(mech, solver, rt, momentum, cst, nx, ny; iters = warmup, gamma, theta_v,
                glen = false, dt_refresh, ncheck, rtol = 1e-8, scale, a1)
    err, iters_used = bp_iterate!(mech, solver, rt, momentum, cst, nx, ny; iters, gamma,
                                  theta_v, glen = true, dt_refresh, ncheck, rtol, scale,
                                  a1, trace, stagnation)
    elapsed = time() - t0
    verticalvelocity!(mech, rt)
    return (; grid, rt, mech, nx, ny, nz, L, err, iters_used, elapsed,
              converged = err <= rtol)
end

#=
## Running Pagos' DIVA solver on the same experiments

`roadmaps/chmy.md` Phase 3. DIVA is depth-integrated exactly like SSA, so — unlike BP — the
library's own per-iteration step, [`pseudo_rate!`](@ref), needs no surgery for periodicity:
the depth-averaged velocity gradients, membrane stress and basal update it assembles are all
either local reads or direct neighbour differences, the same node algebra
[`pseudo_rate!(..., ::MomentumBalance3D, ...)`](@ref) already relies on inside `bp_iterate!`.
What breaks under periodicity is only what [`pseudo_transient!`](@ref) does *outside*
`pseudo_rate!`: the `bc!(…, Neumann())` calls on `velocity.depthaverage_x`/`y` (the same
`Chmy.BoundaryConditions.batch_impl` gap documented above `bp_iterate!`), and the
caller-driven DIVA depth-integrated-viscosity chain ([`diva_update!`](@ref)), whose
`µ̄`/`β_eff` outputs are read one node beyond their own launch range and so need the periodic
seam refreshed explicitly — exactly like BP's viscosity halo.

So [`pseudo_rate!`](@ref) is called unmodified — it is exported precisely so callers can do
this — and only the two things `pseudo_transient!` does around it are replaced: `bc!` becomes
[`periodic_halo!`](@ref), and `diva_update!`'s outputs get an explicit halo refresh before
anything downstream reads them across the seam.

!!! note "Which fields need the seam refreshed, and why membrane stress itself does not"
    `depthaverage_velocitygradients!` writes `x_dy`/`y_dx` (at `ab`) from `depthaverage_x`/`y`
    (periodic already, via the velocity halo below); [`membranestress!`](@ref) then reads
    them at their *own* node, no interpolation, so it needs nothing further. `diva_update!`'s
    `effective_strainrate_diva!`, however, `lerp`s them onto `aa` — a neighbour of an
    already-one-ring-extended field, the same "gradient of a gradient" reach that forces BP's
    terrain correction to halo `x_dz`/`y_dz` — so `depthaverage_x_dy`/`y_dx` need an explicit
    refresh before the *next* `diva_update!` call reads them. `viscosity_depthaveraged` and
    `beta_eff` are `hlerp`/`lerp`ed the same way — in `membranestress!`, in
    `update_basalstress!` and in the Gershgorin bound — so both are haloed right after
    `diva_update!` writes them, before `pseudo_rate!` runs.
=#
function diva_iterate!(mech, solver, rt, momentum, cst, nx, ny; iters, gamma, theta_v,
                       glen::Bool, dt_refresh, ncheck, rtol, scale,
                       trace = nothing, stagnation = true)
    ux, uy = mech.velocity.depthaverage_x, mech.velocity.depthaverage_y
    err, errprev, used = Inf, Inf, 0
    for iter in 1:iters
        used = iter
        if glen
            diva_update!(mech, solver, rt)
            periodic_halo!(mech.material.viscosity_depthaveraged, nx, ny)
            periodic_halo!(mech.friction.beta_eff, nx, ny)
            iter % dt_refresh == 0 && pseudo_dt!(solver, mech, cst, rt)
        end
        pseudo_rate!(mech, cst, rt, momentum, solver; gamma)
        ## Read by the *next* `diva_update!` (see the note above), not by anything in this
        ## iteration — so the refresh can wait until after `pseudo_rate!` has written them.
        periodic_halo!(mech.velocity.depthaverage_x_dy, nx, ny)
        periodic_halo!(mech.velocity.depthaverage_y_dx, nx, ny)
        copyto!(asarray(solver.velocity_x_old), asarray(ux))
        copyto!(asarray(solver.velocity_y_old), asarray(uy))
        pseudo_vel!(asarray(ux), asarray(solver.velocity_x_old),
                    asarray(solver.velocity_x_dt), asarray(solver.dtau_x), theta_v)
        pseudo_vel!(asarray(uy), asarray(solver.velocity_y_old),
                    asarray(solver.velocity_y_dt), asarray(solver.dtau_y), theta_v)
        periodic_halo!(ux, nx, ny); periodic_halo!(uy, nx, ny)
        if iter % ncheck == 0
            err = max(maximum(abs, interior(solver.residual_x)),
                      maximum(abs, interior(solver.residual_y))) / scale
            trace === nothing || push!(trace, (iter, err))
            (isfinite(err) && err > rtol) || break
            ## Same stagnation guard as `bp_iterate!`, for the same frozen-bed reason.
            (stagnation && abs(err - errprev) < 1e-4 * err && iter > 4000) && break
            errprev = err
        end
    end
    return err, used
end

#=
!!! note "The caller fills `friction.beta` (bare), not `beta_eff`"
    Unlike BP and SSA, DIVA *derives* `β_eff` from the bare friction coefficient `β` and the
    viscosity integral `F₂` (Eqs. 19–20 of Robinson et al. 2022) — handing it a pre-filled
    `beta_eff` would let it silently skip the depth-integrated correction ISMIP-HOM is meant
    to exercise. The same numeric value as `BETA_NOSLIP` still gives the no-slip limit: at
    the operating viscosity `F₂ = H/(3µ) ~ 2e-4`, so `β_eff → 1/F₂ ~ 5e3` once `β ≫ 1/F₂`,
    the frozen-bed asymptote (Eq. 20) — bounded by the column's own shear resistance rather
    than by how large `β` is set, which is the physically correct no-slip behaviour.

!!! note "Warm start: a `β_eff` from the guessed `µ`, without calling `diva_update!`"
    The frozen-viscosity warmup stage still needs a `β_eff` to iterate the velocity against,
    but must not run the real Glen continuation on the `u = 0` strain rate — the same
    cold-start trap `run_pagos_bp` avoids by holding `µ` fixed. So the warmup `β_eff` is
    derived by hand from the guessed constant `µ` ([`viscosity_integrals!`](@ref) then
    [`beta_eff_diva!`](@ref), skipping [`update_viscosity!`](@ref) entirely), then held fixed
    exactly like BP's `µ` until the second stage switches `glen` on.
=#
function run_pagos_diva(experiment, L_km; nx = PAGOS_NX[experiment], nz = PAGOS_NZ,
                        beta0 = BETA_NOSLIP, warmup = 2000, iters = 60000, gamma = 0.6,
                        theta_v = 0.6, theta_mu = 0.05, reg = 1e-8, dt_refresh = 20,
                        ncheck = 250, rtol = 1e-5, cfl = 0.9,
                        trace = nothing, stagnation = true)
    L = L_km * 1e3; ω = 2π / L
    ny = experiment === :A ? nx : 4          # exp B has no y dependence
    dx = L / nx
    layering = CorrectedVerticalLayering(T_PAGOS, QuadraticSigmaTransform(T_PAGOS, nz))

    ## Build and fill on the CPU — see "CPU or GPU" above — then move the whole state at once.
    grid_cpu = StaggeredGrid(T_PAGOS, L, ny * dx, dx, dx, layering)
    mech_cpu = MechanicState(grid_cpu)
    bump = experiment === :A ? (x, y) -> sin(ω * x) * sin(ω * y) : (x, y) -> sin(ω * x)
    fill_xy!(mech_cpu.topography.thickness, grid_cpu.grid2d, (x, y) -> 1000 - 500 * bump(x, y))
    fill_xy!(mech_cpu.topography.surface,   grid_cpu.grid2d, (x, y) -> -x * SLOPE)
    fill_xy!(mech_cpu.friction.beta,        grid_cpu.grid2d, (x, y) -> beta0)
    fill_xy!(mech_cpu.material.rate_factor, grid_cpu.grid,   (x, y) -> A_GLEN)
    fill_xy!(mech_cpu.material.viscosity,   grid_cpu.grid,   (x, y) -> MU_GUESS)
    fill_xy!(mech_cpu.material.viscosity_depthaveraged, grid_cpu.grid2d, (x, y) -> MU_GUESS)
    setdata!(mech_cpu.velocity.depthaverage_x, zero(T_PAGOS))
    setdata!(mech_cpu.velocity.depthaverage_y, zero(T_PAGOS))

    grid = USE_GPU ?
        StaggeredGrid(Arch(CUDABackend()), T_PAGOS, L, ny * dx, dx, dx, layering) : grid_cpu
    mech = USE_GPU ? Pagos.Adapt.adapt(CuArray, mech_cpu) : mech_cpu
    rt = Runtime(grid); cst = Constants{T_PAGOS}()
    momentum = DIVAMomentumBalance()
    periodic_halo!(mech.velocity.depthaverage_x, nx, ny)
    periodic_halo!(mech.velocity.depthaverage_y, nx, ny)

    ## `PseudoTransientSolver(grid; …)`, not `(grid, momentum; …)`: DIVA's work arrays live on
    ## `grid.grid2d` like SSA's, not on the column grid BP's constructor builds them on.
    solver = PseudoTransientSolver(grid;
        viscosity_continuation = DIVAViscosityContinuation(T_PAGOS; n_glen = 3, theta_mu,
                                                            strainrate_reg = reg),
        friction_update = ActiveFrictionUpdate(),
        pseudo_timestep = GershgorinPseudoTimeStep(; cfl))

    drivingstress!(mech, cst, rt)
    ## Warmup β_eff from the guessed constant µ — see the note above.
    viscosity_integrals!(mech.material.viscosity_integral_1,
                        mech.material.viscosity_integral_2, mech, rt)
    beta_eff_diva!(mech, rt)
    periodic_halo!(mech.friction.beta_eff, nx, ny)
    ## Same normalizer concept as `run_pagos_bp`'s `scale`, adapted to DIVA's per-*area*
    ## residual (`dotvel!` divides by `ρH`, not `ρ`): the domain-mean thickness stands in for
    ## the per-node `lerp`'d `H` `_driving_rate!` would use, which is more machinery than a
    ## convergence-gate normalizer needs — the ±500 m bed bumps mostly average out of a 1000 m
    ## mean, and `rtol` only has to be within an order of magnitude of "converged" to be useful.
    scale = maximum(abs, interior(mech.stress.driving_x)) /
            (cst.density_ice * mean(interior(mech.topography.thickness)))
    pseudo_dt!(solver, mech, cst, rt)

    t0 = time()
    diva_iterate!(mech, solver, rt, momentum, cst, nx, ny; iters = warmup, gamma, theta_v,
                 glen = false, dt_refresh, ncheck, rtol = 1e-8, scale)
    err, iters_used = diva_iterate!(mech, solver, rt, momentum, cst, nx, ny; iters, gamma,
                                    theta_v, glen = true, dt_refresh, ncheck, rtol, scale,
                                    trace, stagnation)
    elapsed = time() - t0

    ## Post-solve diagnostics: reconstruct the 3D profile ([`velocities3D!`](@ref), Eq. 16) so
    ## surface `vx`/`vz` and bed `Δp` compare against the same quantities BP reports —
    ## `pagos_fields` below reads `mech.velocity.x`/`y`/`x_dx`/`y_dy` and
    ## `mech.material.viscosity` regardless of which momentum balance produced them. The halo
    ## choreography mirrors `bp_iterate!`'s exactly, run once here rather than every iteration.
    velocities3D!(mech, rt, momentum)
    periodic_halo!(mech.velocity.x, nx, ny); periodic_halo!(mech.velocity.y, nx, ny)
    velocitygradients!(mech.velocity, mech.topography.thickness, rt)
    periodic_halo!(mech.velocity.x_dz, nx, ny); periodic_halo!(mech.velocity.y_dz, nx, ny)
    terrain_metric_correction!(mech.velocity, mech.topography.thickness,
                               mech.topography.surface, rt)
    periodic_halo!(mech.velocity.x_dx, nx, ny); periodic_halo!(mech.velocity.y_dy, nx, ny)
    periodic_halo!(mech.velocity.x_dy, nx, ny); periodic_halo!(mech.velocity.y_dx, nx, ny)
    verticalvelocity!(mech, rt)

    return (; grid, rt, mech, nx, ny, nz, L, err, iters_used, elapsed,
              converged = err <= rtol)
end

#=
Diagnostics, mapped onto the reference submissions' own output convention (§ "File format"):
everything is reported at cell centres on the normalized abscissa `x̂ = x/L`. Every field is
pulled to the host with `Array` first — one bulk transfer each, for the same reason the fill
went the other way in bulk.

`Δp` deserves a note, because Blatter-Pattyn does not solve for pressure at all. It does not
have to: the first-order vertical balance is `∂σ_zz/∂z = ρg` with `σ_zz = 0` at the surface,
so `σ_zz = -ρg(s-z)`, and with `σ_zz = -p_I + τ_zz` that gives `p_I = p_H + τ_zz` exactly.
Hence

```math
\Delta p = p_I - p_H = \tau_{zz} = 2\mu\dot\varepsilon_{zz} = -2\mu\,(u_x + v_y),
```

evaluated in the bottom layer — every factor of which BP already carries. `τ_xz(z_b)` is the
bed flux the solver imposes, `β u_b`, i.e. `stress.base_x`.
=#
function pagos_fields(r)
    (; mech, nx, ny, nz, L) = r
    host(f) = Array(interior(f))
    vxf = host(mech.velocity.x); vyf = host(mech.velocity.y)
    w   = host(mech.velocity.z); μ = host(mech.material.viscosity)
    ux  = host(mech.velocity.x_dx); vy = host(mech.velocity.y_dy)
    bx  = host(mech.stress.base_x); by = host(mech.stress.base_y)
    dx  = L / nx
    ctr(i) = (i - 0.5) * dx / L
    x = [ctr(i) for i in 1:nx]; y = [ctr(j) for j in 1:ny]
    ## acx/acy → aa by averaging the two faces; w is already horizontally at aa.
    vxs = [(vxf[i, j, nz] + vxf[i + 1, j, nz]) / 2 for i in 1:nx, j in 1:ny]
    vys = [(vyf[i, j, nz] + vyf[i, j + 1, nz]) / 2 for i in 1:nx, j in 1:ny]
    vzs = [w[i, j, nz + 1] for i in 1:nx, j in 1:ny]
    txz = [(bx[i, j, 1] + bx[i + 1, j, 1]) / 2000 for i in 1:nx, j in 1:ny]   # kPa
    tyz = [(by[i, j, 1] + by[i, j + 1, 1]) / 2000 for i in 1:nx, j in 1:ny]
    dp  = [-2 * μ[i, j, 1] * (ux[i, j, 1] + vy[i, j, 1]) / 1000 for i in 1:nx, j in 1:ny]
    return (; x, y, vx = vxs, vy = vys, vz = vzs, txz, tyz, dp)
end

#=
## Plot style

Two things to encode. `L` is an **ordered** variable, not a set of unrelated categories, so
the five domain lengths are coloured by a sequential ramp — the reader should see "shorter
domain" as a direction, not look the colour up each time. The ramp is truncated short of
`plasma`'s pale yellow end, which would be unreadable as a 2px line on white.

Model identity takes the *fill vs. line* channel: the full-Stokes ensemble is a translucent
band, Pagos BP is a solid line on top of it. Reading "is BP inside the band?" is then a
single visual question at every `x̂`, which is the question the whole comparison is for.
=#
set_theme!(theme_latexfonts())
update_theme!(Axis = (xgridcolor = (:black, 0.06), ygridcolor = (:black, 0.06),
                      xgridwidth = 1, ygridwidth = 1,
                      topspinevisible = false, rightspinevisible = false,
                      spinewidth = 0.8, titlesize = 15))

const LCOLORS = [cgrad(:plasma)[t] for t in range(0.05, 0.80; length = length(LENGTHS_KM))]
const LLABELS = string.(LENGTHS_KM)
const XQ = collect(range(0, 1; length = 201))

#=
Periodic linear interpolation onto the common abscissa the band needs. The submissions'
abscissae differ (cell centres on `(0, 1)` against nodes on `[0, 1]`), and the field is
periodic, so the table is extended by one wrapped point at each end rather than clamped —
clamping would flatten the profile over the first and last half-cell, exactly where the bed
crest sits for some `L`. A submission whose grid already contains both endpoints produces a
zero-width interval at the seam, which is why the zero denominator is guarded rather than
assumed impossible.
=#
function interp_periodic(x, v, xq = XQ)
    xs = vcat(x[end] - 1, x, x[1] + 1)
    vs = vcat(v[end], v, v[1])
    out = similar(xq, Float64)
    for (k, q) in enumerate(xq)
        i = clamp(searchsortedlast(xs, q), 1, length(xs) - 1)
        d = xs[i + 1] - xs[i]
        out[k] = d == 0 ? vs[i] : (1 - (q - xs[i]) / d) * vs[i] + ((q - xs[i]) / d) * vs[i + 1]
    end
    return out
end

#=
The `Δp` sign convention, as a single signed number per submission: `Δp` varies along the
profile essentially as the bed's fundamental mode, so projecting it onto `cos(2πx̂)` isolates
the convention and discards the magnitude. Grid-independent by construction (it is an average
over whatever abscissa the submission used), which is what lets it work across the whole
ensemble rather than needing a hard-coded flip per model.
=#
dp_orientation(p) = sign(mean(p.dp .* cos.(2π .* p.x)))

#=
### Not every submission's column 8 is the same quantity

Sign is not the only thing that varies in `Δp`. Across the five full-Stokes submissions the
velocity columns agree to 0.5–2.6 % at every domain length, but the `Δp` *amplitudes* span a
factor of 20: at `L = 160 km`, `rhi3`/`oga1`/`rhi1` all report ±3.3 kPa while `aas2` reports
±60. The means are ~0 in every case, so this is not a constant offset — it is
a different quantity under the same column heading.

Which three are right is not a matter of taste, because `Δp` has a limit to satisfy: it is
**identically zero under the shallow-ice approximation**, so it must collapse as `L` grows and
the geometry becomes shallow-ice-like. Measured as the ratio of the `L = 160 km` amplitude to
the `L = 10 km` one:

| | rhi3 | oga1 | aas2 | rhi1 | (cma1) |
|---|---|---|---|---|---|
| amp(160)/amp(10) | 0.089 | 0.088 | **0.90** | 0.089 | **0.95** |

Three collapse by a factor of 11; `aas2` is essentially flat in `L`, which no pressure
anomaly can be. (`cma1`, dropped from `MODELS` for resolution, failed the same way.)

So the test is the physics, not the peer group — a submission whose `Δp` does not
decay towards the shallow-ice limit is reporting something else (plausibly a stress component
or an unsubtracted overburden term), and averaging it into a min–max band would inflate the
band 40-fold and drown the quantity being measured.

Applied to the `Δp` band and the peak-`|Δp|` summary only; `vx`, `vz` and `τxz` use the whole
ensemble, where all five agree. Excluded members are named in the printed table, not dropped
quietly.

The test is per-experiment, and that turned out to matter when experiment B was still being
read: **all five submissions passed it there**, reporting `Δp = 3.3–4.2 kPa` at
`L = 160 km`, `aas2` included. Whatever goes wrong in their experiment-A files does not go
wrong in their B files — which is the strongest evidence that it is a reporting artifact of
the 2D output rather than a property of the model. Nothing here needs to decide which; the
criterion is applied where it fires and nowhere else.
=#
function dp_ensemble_mask(refs; ratio = 0.3)
    amp(p) = (maximum(p.dp) - minimum(p.dp)) / 2
    long, short = refs[1], refs[end]      # LENGTHS_KM runs longest → shortest
    return [amp(long[m]) < ratio * amp(short[m]) for m in eachindex(long)]
end

#=
The full-Stokes envelope at one domain length: min and max across the ensemble, on `XQ`.
`Δp` is the one quantity whose members are sign-normalised first (see the warning at the top);
`flip` carries the per-model sign so that decision stays visible at the call site instead of
being buried here.
=#
function fs_envelope(profiles, key; flip = ones(length(profiles)))
    M = reduce(hcat, (interp_periodic(p.x, getproperty(p, key)) .* s
                      for (p, s) in zip(profiles, flip)))
    return vec(minimum(M; dims = 2)), vec(maximum(M; dims = 2))
end

const PANEL_KEYS = ((:vx,  "surface velocity vx(zs)  (m/yr)"),
                    (:vz,  "surface velocity vz(zs)  (m/yr)"),
                    (:txz, "basal shear stress τxz(zb)  (kPa)"),
                    (:dp,  "pressure anomaly Δp(zb)  (kPa)"))

#=
Four panels rather than one with four y-scales: the quantities differ by three orders of
magnitude and share only their abscissa, so a twin axis would be unreadable *and* would
invite false slope comparisons.

`refs[i]` is the vector of reference profiles at `LENGTHS_KM[i]`; `pagos[i]` is the single
Pagos profile there, or `nothing` to plot the reference ensemble alone. `pagos2` is a second
Pagos curve (DIVA, run alongside BP) plotted **dashed** rather than with a second colour
channel — colour is already spent on domain length, so a second model reuses it and encodes
its identity in linestyle instead, exactly the fill-vs-line split the module docs above
describe for the reference band vs. `pagos`.

`FS_ALPHA` is shared with the Legend swatches below so the two never drift apart; it sits
above the fully-transparent end of the range precisely because five overlapping bands at very
low alpha wash out to indistinguishable pale colour — this value is picked to keep the
per-length band readable while still letting the Pagos lines on top show through.
=#
const FS_ALPHA = 0.38

function profile_figure(refs, pagos, suptitle, nmodels;
                        flips = fill(ones(nmodels), length(LENGTHS_KM)),
                        dpmask = trues(nmodels), pagos2 = nothing,
                        pagos_label = "Pagos BP", pagos2_label = "Pagos DIVA")
    fig = Figure(size = (1250, 820))
    Label(fig[0, 1:3], suptitle; fontsize = 19, font = :bold, padding = (0, 0, 6, 0))
    for (p, (key, lab)) in enumerate(PANEL_KEYS)
        row, col = fldmod1(p, 2)
        ## The `Δp` panel states both of its departures from the other three — the band is
        ## sign-normalised and restricted to the submissions whose `Δp` satisfies the
        ## shallow-ice limit. Unlabelled, either would read as a physics spread.
        note = key === :dp ?
            "band: $(count(dpmask))/$nmodels submissions (shallow-ice limit), " *
            "sign-normalised; Pagos as computed" : ""
        ax = Axis(fig[row, col]; xlabel = "normalized x  (x / L)", ylabel = lab,
                  title = note, titlesize = 12, titlecolor = :firebrick)
        for (i, _) in enumerate(LENGTHS_KM)
            members = key === :dp ? refs[i][dpmask] : refs[i]
            flip = key === :dp ? flips[i][dpmask] : ones(length(members))
            lo, hi = fs_envelope(members, key; flip)
            band!(ax, XQ, lo, hi; color = (LCOLORS[i], FS_ALPHA))
        end
        if pagos !== nothing
            for (i, _) in enumerate(LENGTHS_KM)
                lines!(ax, pagos[i].x, getproperty(pagos[i], key);
                       color = LCOLORS[i], linewidth = 2)
            end
        end
        if pagos2 !== nothing
            for (i, _) in enumerate(LENGTHS_KM)
                lines!(ax, pagos2[i].x, getproperty(pagos2[i], key);
                       color = LCOLORS[i], linewidth = 2, linestyle = :dash)
            end
        end
        xlims!(ax, 0, 1)
    end
    model_elements = [PolyElement(color = (:gray30, FS_ALPHA)),
                      LineElement(color = :gray30, linewidth = 2)]
    model_labels = ["full-Stokes spread\n($nmodels submissions)", pagos_label]
    if pagos2 !== nothing
        push!(model_elements, LineElement(color = :gray30, linewidth = 2, linestyle = :dash))
        push!(model_labels, pagos2_label)
    end
    Legend(fig[1:2, 3],
           [[PolyElement(color = (c, FS_ALPHA)) for c in LCOLORS], model_elements],
           [LLABELS, model_labels],
           ["domain length\nL (km)", "model"]; framevisible = false)
    colsize!(fig.layout, 3, Relative(0.15))
    return fig
end

#=
How the solution scales with `L`, summarising in three numbers what the profile panels show
one at a time, on a log abscissa (the domain lengths are a geometric sequence). Three separate
small axes rather than one with three bands: the quantities have unrelated units and ranges.

 - **mean `vx`** — how fast the slab flows on average;
 - **profile amplitude** `max − min` — how strongly the bed still prints through to the
   surface, which is the quantity longitudinal stresses destroy as `L` shrinks;
 - **peak `|Δp|`** — the departure of bed pressure from hydrostatic, identically zero under
   the shallow-ice approximation and therefore the cleanest scalar measure of how much
   higher-order physics the experiment is exercising.
=#
const SUMMARIES = (("mean vx(zs)  (m/yr)", "profile mean", p -> mean(p.vx)),
                   ("max − min of vx(zs)  (m/yr)", "bed signature at the surface",
                    p -> maximum(p.vx) - minimum(p.vx)),
                   ("max |Δp(zb)|  (kPa)", "departure from hydrostatic (0 under SIA)",
                    p -> maximum(abs, p.dp)))

function scaling_figure(refs, pagos, suptitle; dpmask = trues(length(refs[1])))
    fig = Figure(size = (1350, 400))
    Label(fig[0, 1:3], suptitle; fontsize = 19, font = :bold, padding = (0, 0, 6, 0))
    Lf = Float64.(LENGTHS_KM)
    for (s, (ylab, subtitle, f)) in enumerate(SUMMARIES)
        ## The third summary is peak |Δp|, so it takes the restricted ensemble — same reason
        ## the `Δp` panel does.
        members(i) = s == 3 ? refs[i][dpmask] : refs[i]
        ax = Axis(fig[1, s]; xscale = log10, xlabel = "domain length L (km)", ylabel = ylab,
                  title = s == 3 ? subtitle * "; $(count(dpmask)) submissions" : subtitle,
                  titlesize = 13, xticks = (Lf, string.(LENGTHS_KM)))
        lo = [minimum(f(p) for p in members(i)) for i in eachindex(Lf)]
        hi = [maximum(f(p) for p in members(i)) for i in eachindex(Lf)]
        band!(ax, Lf, lo, hi; color = (:gray30, 0.25))
        if pagos !== nothing
            v = [f(pagos[i]) for i in eachindex(Lf)]
            lines!(ax, Lf, v; color = LCOLORS[end], linewidth = 2)
            scatter!(ax, Lf, v; color = LCOLORS[end], markersize = 8)
        end
    end
    return fig
end

#=
The table the figures are read against. The **full-Stokes spread** column is the one to read
first: the ensemble max−min in mean `vx`, as a percentage of the ensemble mean. That is the
floor on how finely any approximate model can be judged at each domain length — a Pagos
result inside it is indistinguishable from full Stokes given this reference set, and a
disagreement below it is not evidence of anything. `Pagos−FS` is measured against the
ensemble *mean*, and is only meaningful when it exceeds the spread.
=#
function results_table(experiment, refs, pagos, labels; dpmask = trues(length(labels)))
    println()
    println("experiment $experiment — mean / max surface vx (m/yr), peak |Δp| (kPa)")
    print(rpad("L (km)", 8))
    for lab in labels
        print("| ", rpad(lab, 23))
    end
    print("| ", rpad("Pagos BP", 23), "|  FS spread |  Pagos-FS")
    println()
    println("-"^(8 + 25 * (length(labels) + 1) + 24))
    for (i, L) in enumerate(LENGTHS_KM)
        print(rpad(L, 8))
        for (m, p) in enumerate(refs[i])
            ## `†` = this submission's `Δp` fails the shallow-ice-limit test and is excluded
            ## from the band; its velocity columns are still in the ensemble.
            print(@sprintf("| %7.2f %7.2f %6.1f%s", mean(p.vx), maximum(p.vx),
                           maximum(abs, p.dp), dpmask[m] ? " " : "†"))
        end
        means = [mean(p.vx) for p in refs[i]]
        fsmean = mean(means)
        spread = 100 * (maximum(means) - minimum(means)) / fsmean
        if pagos === nothing
            print(@sprintf("| %-23s", "—"))
            println(@sprintf("|  %8.2f%% |        —", spread))
        else
            p = pagos[i]
            flag = hasproperty(p, :converged) && !p.converged ? "!" : " "
            print(@sprintf("| %7.2f %7.2f %6.1f%s", mean(p.vx), maximum(p.vx),
                           maximum(abs, p.dp), flag))
            println(@sprintf("|  %8.2f%% |  %+7.2f%%", spread,
                             100 * (mean(p.vx) - fsmean) / fsmean))
        end
    end
    println("\n! = solve did not reach the residual tolerance; see the Pagos section above.")
    all(dpmask) ||
        println("† = Δp excluded from the band (fails the shallow-ice limit): " *
                join(labels[.!dpmask], ", "))
    println("models: $(join(labels, ", ")) — ISMIP-HOM (Pattyn et al., 2008)")
end

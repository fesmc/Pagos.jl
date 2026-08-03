#=

# [Antarctic momentum-balance comparisons](@id ais_momentum)

The merged successor to the earlier `ais-pt.jl` (SSA only) and `ais-diva.jl` (DIVA only,
never run — no data was available when it was written). One real 761×761, 8 km Antarctic
Ice Sheet geometry, one Yelmo restart file supplying *both* the 2D fields the SSA limit
needs and the genuinely 3D viscosity/velocity DIVA needs, and five comparisons built on top
of the same converged state:

 1. SSA vs. DIVA, both run in Pagos, on identical geometry/friction/viscosity.
 2. Pagos DIVA vs. Yelmo's own DIVA solve (surface speed only, RMSE — not the full 3D
    profile: Yelmo's own vertical grid does not need to match Pagos' for this).
 3. `FixedTuning()` (untried defaults) vs. [`AutotunedDynamicRelaxation`](@ref), for DIVA —
    the SSA version of this comparison is already the headline result of
    `roadmaps/PT-autotune.md` Phase 2 (6.3× on this exact geometry); this is its DIVA
    counterpart, which nobody has run before.
 4. CPU vs. GPU wall-clock, same solve.
 5. Float64 vs. Float32 wall-clock and accuracy.

None of the five is a validated simulation — see `ais-pt.jl`'s own framing, which still
applies: this is a code-state check on real geometry, not a benchmark paper.

!!! warning "Comparison 1 is not quite apples-to-apples on friction"
    The SSA run uses Yelmo's own `beta_eff` directly (as [`ais-pt.jl`](@ref ais_pt) always
    did) — but that field is *already* Yelmo's own DIVA vertical-shear correction, not a
    bare SSA friction coefficient. The DIVA run instead starts from the bare `beta` and
    derives its own `β_eff` via [`diva_update!`](@ref), from Pagos' own 3D viscosity and
    geometry — which will not exactly reproduce Yelmo's `beta_eff`, not least because of the
    layering mismatch noted below. So comparison 1 is "SSA fed Yelmo's DIVA-corrected
    friction" vs. "DIVA deriving its own correction from the same raw inputs" — informative
    about whether Pagos' DIVA path does something structurally different from SSA at all,
    but not a controlled ablation of the friction law alone.

!!! note "Vertical layering: index-matched to Yelmo, not value-matched"
    Yelmo's `zeta`/`zeta_ac` do not correspond to [`QuadraticSigmaTransform`](@ref)'s own
    derivation (checked numerically — different values at every interior level), and
    reproducing them exactly would mean bypassing the parametric transform for a
    `CorrectedVerticalLayering` built straight from Yelmo's arrays. Not worth it here: this
    script uses a plain `QuadraticSigmaTransform` grid with Yelmo's own layer *count*
    (`nz = 11`) and copies Yelmo's per-layer viscosity in by index, not by matching sigma
    value. That is why comparison 2 below is surface-only (a quantity insensitive to the
    exact interior discretization) rather than a profile comparison.

!!! warning "Memory budget: written for a 16 GB machine, not a cluster node"
    A single [`MechanicState`](@ref) on this grid (761×761×11) is on the order of **2 GB**
    — `StressState`, `VelocityState` and `StrainRateState` between them hold ~35 column
    (`AA3`/`AB3`/…) fields at `(761, 761, 11)` Float64, ~51 MB each. The first version of
    this script kept *two* of them (`mech_ssa`, `mech_diva`) alive simultaneously, plus every
    raw NetCDF array at Float64 (including the full 3D `ux`/`uy`, never actually needed at
    more than their surface layer) — comfortably enough to exhaust 16 GB and get the whole
    process killed with no Julia-level error at all (confirmed: the run that triggered this
    rewrite left no completion record, and `dmesg`/`free -h` are consistent with an OOM
    kill, not a Pagos bug). Three changes fix this, applied throughout below:

     1. **Every solve lives inside [`run_solve`](@ref)**, which returns only a small 2D
        speed field and a handful of scalars — never the `MechanicState` itself. `mech` is
        local to that function and becomes collectible the moment it returns; an explicit
        `GC.gc()` right after each call forces that to happen *before* the next multi-GB
        allocation, rather than trusting Julia's lazy GC heuristic to keep up.
     2. **Loaded fields use the narrowest safe type**, not a blanket cast. Checked against
        the file's own value ranges first (`visc` peaks at `9.6e9`, `beta` at `1.3e9` — both
        would silently become `Inf` under Float16, whose max magnitude is `65504`): bounded
        2D fields (`H_ice`, `z_srf`, `z_bed`, `f_grnd`, Yelmo's surface speed) are safe at
        Float16; viscosity and friction, which are not, stay Float32 — still half the
        memory of the Float64 this script started with.
     3. **`ux`/`uy` are read as a 2D surface slice directly from the file**
        (`ds["ux"][:, :, end, 1]`), not the full `(761, 761, 11)` column — comparison 2 only
        ever needed the surface, so the other 10 layers are never brought into memory at
        all, and the read happens where that comparison is built, not upfront.

    `ram()` below prints free/total system memory (not just Julia's own heap) at every
    section boundary, so the effect of each of these is visible while the script runs
    rather than asserted.

=#

using Pagos, NCDatasets, CairoMakie, Statistics, Printf

restart_file = "/home/jan/pCloudSync/PhD/Projects/Ice-Sheet-Modelling/ice-data-pagos/yelmo_restart_ais_8km.nc"

ram() = @printf("  free RAM: %s / %s\n",
               Base.format_bytes(Sys.free_memory()), Base.format_bytes(Sys.total_memory()))

println("Before loading anything:"); ram()

#=
## Load geometry, at the narrowest type each field's own value range allows

`load2d`/`load3d` read one variable, drop the trailing size-1 `time` dimension, and cast
directly to `T` in a single pass — no Float64 intermediate that then gets thrown away.
Value ranges checked against the file beforehand (see the memory note above): `H_ice`,
`z_srf`, `z_bed`, `f_grnd` and Yelmo's own depth-averaged/surface speeds all stay
comfortably under Float16's `65504` magnitude cap; `visc_bar`/`visc`/`beta`/`beta_eff` do
not (up to `9.6e9`), so those load as Float32.
=#

load2d(ds, name, T) = T.(dropdims(ds[name][:, :, :]; dims = 3))
load3d(ds, name, T) = T.(dropdims(ds[name][:, :, :, :]; dims = 4))

xc, yc, H_ice, z_srf, z_bed, f_grnd, visc_bar, beta, beta_eff, visc3d, zeta =
NCDataset(restart_file) do ds
    (Float64.(ds["xc"][:]), Float64.(ds["yc"][:]),
     load2d(ds, "H_ice", Float16), load2d(ds, "z_srf", Float16), load2d(ds, "z_bed", Float16),
     load2d(ds, "f_grnd", Float16),
     load2d(ds, "visc_bar", Float32), load2d(ds, "beta", Float32), load2d(ds, "beta_eff", Float32),
     load3d(ds, "visc", Float32),
     Float64.(ds["zeta"][:]))
end

nx, ny = length(xc), length(yc)
nz = length(zeta)
dx = (xc[2] - xc[1]) * 1e3
dy = (yc[2] - yc[1]) * 1e3
lx, ly = nx * dx, ny * dy

println("After loading NetCDF fields (Float16/Float32, ux/uy surface deferred to §2):")
ram()

#=
## Build the grid, topography state, and ice mask

One [`StaggeredGrid`](@ref) with a real column (`nz = 11`, [`QuadraticSigmaTransform`](@ref)
— see the layering note above), shared by every comparison below: [`TopographicState`](@ref)
and the mask are built once, since nothing about them depends on which momentum balance or
element type a given solve uses. The mask follows [`ais-pt.jl`](@ref ais_pt)'s own reasoning
exactly — [`momentum_mask!`](@ref)'s `is_momentum_solved`, not `is_ice`, to exclude the ~48
detached iceberg cells a force balance cannot be posed on.
=#

layering = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, nz))
grid = StaggeredGrid(Float64, lx, ly, dx, dy, layering)
rt   = Runtime(grid)
topo = TopographicState(grid)
cst  = Constants{Float64}()

function fill_from_grid!(f, data)
    ni, nj = size(data)
    for k in axes(interior(f), 3), j in -1:(nj + 2), i in -1:(ni + 2)
        ic, jc = clamp(i, 1, ni), clamp(j, 1, nj)
        f[i, j, k] = data[ic, jc]
    end
    return f
end

# As `fill_from_grid!`, but `data` is genuinely 3D (`nx, ny, nz`) and each layer `k` of `f`
# reads the matching layer of `data` — not a 2D field broadcast down every layer.
function fill_from_grid3d!(f, data)
    ni, nj, nk = size(data)
    for k in 1:nk, j in -1:(nj + 2), i in -1:(ni + 2)
        ic, jc = clamp(i, 1, ni), clamp(j, 1, nj)
        f[i, j, k] = data[ic, jc, k]
    end
    return f
end

fill_from_grid!(topo.thickness.ice, H_ice)
fill_from_grid!(topo.mask.is_grounded, f_grnd .> 0)
icemasks!(topo, rt)
momentum_mask!(topo, rt)
mask = IceMask(topo.mask.is_momentum_solved)

n_detached = count(asarray(topo.mask.is_ice) .& .!asarray(topo.mask.is_momentum_solved))
@show n_detached

#=
## Unit conversion, in place

`Constants` is SI throughout; Yelmo's viscosity/friction are Pa yr / Pa yr m⁻¹.
`@.`-fused broadcasts write back into the same array (`ifelse` reads every element before
any of them are overwritten, so this is safe despite LHS and RHS naming the same binding),
so converting and off-ice-flooring these never allocates a second full-size copy — the
narrow-type load above would otherwise be undone immediately by a same-sized Float64-scale
transient.
=#

spy = cst.seconds_per_year
@. visc_bar *= spy
@. visc3d   *= spy
@. beta     *= spy
@. beta_eff *= spy
visc_off_ice   = maximum(visc_bar)
visc3d_off_ice = maximum(visc3d)
@. visc_bar = ifelse(H_ice > 0, visc_bar, visc_off_ice)
@. visc3d   = ifelse(H_ice > 0, visc3d, visc3d_off_ice)

on_ice(a) = ifelse.(H_ice .> 0, a, NaN)

println("After unit conversion (in place, no new full-size arrays):"); ram()

#=
## `run_solve`: the one place every comparison below drives a solve

Builds a [`MechanicState`](@ref) **locally** — it is never returned, so it becomes eligible
for garbage collection the instant this function returns, rather than accumulating at the
caller's top level the way `mech_ssa`/`mech_diva` did in the first version of this script.
Only a small 2D speed field (`nx × ny`, ~2 MB at Float32, not the ~2 GB `mech` it came from)
and a few scalars survive.

`T` is **derived from `grid` itself** (`eltype(grid.grid)`), not taken as a separate
argument: `MechanicState(grid)` would silently keep allocating Float64 fields regardless of
any `T` passed in independently, since it reads the element type off `grid`, not off a
caller's say-so — a real trap for §5's Float64-vs-Float32 comparison, caught before it was
ever run. §5 gets a genuine Float32 grid by building one, not by passing `T = Float32` to a
Float64 grid. §4 will do the analogous thing for CPU vs. GPU via `Arch`.
=#

function run_solve(momentum, grid, rt, mask; solver_kwargs...)
    T = eltype(grid.grid)
    mech = MechanicState(grid)
    cst_T = Constants{T}()

    fill_from_grid!(mech.topography.thickness, T.(H_ice))
    fill_from_grid!(mech.topography.surface, T.(z_srf))
    fill_from_grid!(mech.material.viscosity_depthaveraged, T.(visc_bar))
    fill_from_grid!(mech.friction.beta_eff, T.(beta_eff))
    fill_from_grid!(mech.friction.beta, T.(beta))
    fill_from_grid3d!(mech.material.viscosity, T.(visc3d))
    setdata!(mech.velocity.depthaverage_x, zero(T))
    setdata!(mech.velocity.depthaverage_y, zero(T))

    solver = PseudoTransientSolver(grid; solver_kwargs...)
    momentum isa DIVAMomentumBalance && diva_update!(mech, solver, rt, mask)

    elapsed = @elapsed result = pseudo_transient!(mech, cst_T, solver, rt, momentum, mask)

    speed = Float32.(spy .* sqrt.(
        (@views (interior(mech.velocity.depthaverage_x)[1:(end - 1), :, 1] .+
                 interior(mech.velocity.depthaverage_x)[2:end, :, 1]) ./ 2) .^ 2 .+
        (@views (interior(mech.velocity.depthaverage_y)[:, 1:(end - 1), 1] .+
                 interior(mech.velocity.depthaverage_y)[:, 2:end, 1]) ./ 2) .^ 2))

    mech = nothing
    GC.gc()
    return (; speed, elapsed, iterations = result.iterations, converged = result.converged,
           residual = result.residual)
end

const SOLVER_KWARGS = (
    abstol = 1e-3, maxiter = 2000, ncheck = 20, printout_every = 50,
    pseudo_timestep = GershgorinPseudoTimeStep(cfl = 0.99),
    convergence = ScaledResidual(),
    friction_update = ActiveFrictionUpdate(),
    tuning = AutotunedDynamicRelaxation(),
)

#=
## 1. SSA vs. DIVA, both in Pagos

Same geometry, same friction, same viscosity magnitude — the only difference is whether the
momentum balance resolves vertical shear. Each call to [`run_solve`](@ref) leaves at most
one `MechanicState` alive at a time.
=#

ssa  = run_solve(SSAMomentumBalance(), grid, rt, mask; SOLVER_KWARGS...)
println("SSA:  ", (; ssa.converged, ssa.iterations, ssa.elapsed, ssa.residual)); ram()

diva = run_solve(DIVAMomentumBalance(), grid, rt, mask; SOLVER_KWARGS...)
println("DIVA: ", (; diva.converged, diva.iterations, diva.elapsed, diva.residual)); ram()

#=
A shared 0–900 m/yr colour scale can hide a real difference between two panels that "look
the same" — the two-panel comparison alone found SSA and DIVA visually near-identical here,
which needed a number, not an impression, before trusting it. `on_ice_mask` excludes ice-free
cells from the statistics (`on_ice` alone just turns them `NaN`, which `filter(!isnan, ...)`
below still has to strip).
=#

on_ice_mask = H_ice .> 0
diff = on_ice(diva.speed .- ssa.speed)
diff_vals = filter(!isnan, diff)
rmse = sqrt(mean(abs2, diff_vals))
rel  = diff_vals ./ max.(filter(!isnan, on_ice(ssa.speed)), 1.0)   # avoid /0 on stagnant ice
# `@printf`'s format string must be a single literal — it is parsed into the generated
# call at macro-expansion time, before a `*`-concatenated string could ever be evaluated.
# `median`, unlike `mean`, has no `f, itr` two-argument method — needs the mapped
# collection built first.
@printf("SSA vs DIVA (on-ice, %d cells): RMSE = %.4g m/yr, mean|Δ| = %.4g m/yr, median rel. diff = %.4g%%, max|Δ| = %.4g m/yr\n",
       length(diff_vals), rmse, mean(abs, diff_vals), 100 * median(abs.(rel)), maximum(abs, diff_vals))

set_theme!(theme_latexfonts())
fig1 = Figure(size = (1650, 620))
crange1 = (0, quantile(filter(!isnan, on_ice(diva.speed)), 0.995))
for (col, (data, title)) in enumerate(((on_ice(ssa.speed), "Pagos SSA"),
                                       (on_ice(diva.speed), "Pagos DIVA")))
    ax = Axis(fig1[1, col], xlabel = "x (km)", ylabel = col == 1 ? "y (km)" : "",
        aspect = DataAspect(), title = title)
    hm = heatmap!(ax, xc, yc, data; colorrange = crange1)
    col == 2 && Colorbar(fig1[1, 3], hm, label = "speed (m/yr)")
end
drange = maximum(abs, diff_vals)
ax3 = Axis(fig1[1, 4], xlabel = "x (km)", aspect = DataAspect(), title = "DIVA - SSA")
hm3 = heatmap!(ax3, xc, yc, diff; colorrange = (-drange, drange), colormap = :RdBu)
Colorbar(fig1[1, 5], hm3, label = "Δ speed (m/yr)")
save("ais-momentum-ssa-vs-diva.png", fig1)
fig1

#=
## 5. Float64 vs. Float32

A genuinely Float32 [`StaggeredGrid`](@ref) — passing `T = Float32` to [`run_solve`](@ref)
alone would *not* do this (see its docstring note); `MechanicState`'s field types come from
the grid, so a separate grid is what actually changes the arithmetic. `_sigma_axis` casts
`layering`'s values to whichever `T` the grid asks for, so the same `Float64`-built
`layering` object is reused rather than needing its own Float32 copy. Reuses the Float64
DIVA result already computed in §1 as the reference, rather than re-solving it.
=#

grid32 = StaggeredGrid(Float32, lx, ly, dx, dy, layering)
rt32   = Runtime(grid32)

diva32 = run_solve(DIVAMomentumBalance(), grid32, rt32, mask; SOLVER_KWARGS...)
println("DIVA Float32: ", (; diva32.converged, diva32.iterations, diva32.elapsed, diva32.residual))
ram()

diff32 = on_ice(Float32.(diva.speed) .- diva32.speed)
diff32_vals = filter(!isnan, diff32)
@printf("Float64 vs Float32 DIVA (on-ice): RMSE = %.4g m/yr, max|Δ| = %.4g m/yr, elapsed %.3gs vs %.3gs (%.2gx)\n",
       sqrt(mean(abs2, diff32_vals)), maximum(abs, diff32_vals), diva.elapsed, diva32.elapsed,
       diva.elapsed / diva32.elapsed)
println("  (elapsed timings run in one process: Float64 kernels were already JIT-compiled ",
       "by §1's SSA/DIVA runs above, Float32 ones compile here for the first time — some ",
       "of the Float32 wall-clock is compilation, not arithmetic. Not a clean throughput ",
       "comparison as measured; a fair one needs each precision timed in its own fresh ",
       "process, or a discarded warm-up solve before the timed one.")

fig2 = Figure(size = (1150, 620))
crange2 = (0, quantile(filter(!isnan, on_ice(diva.speed)), 0.995))
for (col, (data, title)) in enumerate(((on_ice(diva.speed), "DIVA, Float64"),
                                       (on_ice(diva32.speed), "DIVA, Float32")))
    ax = Axis(fig2[1, col], xlabel = "x (km)", ylabel = col == 1 ? "y (km)" : "",
        aspect = DataAspect(), title = title)
    hm = heatmap!(ax, xc, yc, data; colorrange = crange2)
    col == 2 && Colorbar(fig2[1, 3], hm, label = "speed (m/yr)")
end
save("ais-momentum-f64-vs-f32.png", fig2)
fig2

#=
## 3. `FixedTuning()` defaults vs. `AutotunedDynamicRelaxation()`, for DIVA

The DIVA counterpart of `roadmaps/PT-autotune.md` Phase 2's headline SSA result (6.3× on
this exact geometry, `gamma = 0.2` found by a hand scan there). No DIVA hand-tuned value
exists anywhere to reuse, and a genuine hand scan means several full solves before the
comparison even starts (per the earlier discussion) — so this compares the **untuned
library default** `FixedTuning()` (`theta_v = 0.6, gamma = 1`) against the autotuner,
rather than simulating a hand search. `diva` (§1) already *is* the
`AutotunedDynamicRelaxation` reference — no need to solve it twice.
=#

diva_fixed = run_solve(DIVAMomentumBalance(), grid, rt, mask;
                       (; SOLVER_KWARGS..., tuning = FixedTuning(), maxiter = 5000)...)
println("DIVA, FixedTuning() default: ",
       (; diva_fixed.converged, diva_fixed.iterations, diva_fixed.elapsed, diva_fixed.residual))
ram()

@printf("FixedTuning() default vs AutotunedDynamicRelaxation(), DIVA: %d vs %d iterations, %.3gs vs %.3gs (%.2gx)\n",
       diva_fixed.iterations, diva.iterations, diva_fixed.elapsed, diva.elapsed,
       diva_fixed.elapsed / diva.elapsed)

fig3 = Figure(size = (1150, 620))
crange3 = (0, quantile(filter(!isnan, on_ice(diva.speed)), 0.995))
for (col, (data, title)) in enumerate(((on_ice(diva_fixed.speed), "DIVA, FixedTuning()"),
                                       (on_ice(diva.speed), "DIVA, AutotunedDynamicRelaxation()")))
    ax = Axis(fig3[1, col], xlabel = "x (km)", ylabel = col == 1 ? "y (km)" : "",
        aspect = DataAspect(), title = title)
    hm = heatmap!(ax, xc, yc, data; colorrange = crange3)
    col == 2 && Colorbar(fig3[1, 3], hm, label = "speed (m/yr)")
end
Label(fig3[2, 1:3],
    "FixedTuning(): $(diva_fixed.iterations) iter, $(round(diva_fixed.elapsed, digits = 2))s, " *
    "converged = $(diva_fixed.converged), residual = $(round(diva_fixed.residual, sigdigits = 3))" *
    "   |   Autotuned: $(diva.iterations) iter, $(round(diva.elapsed, digits = 2))s, " *
    "converged = $(diva.converged), residual = $(round(diva.residual, sigdigits = 3))" *
    "   |   speedup = $(round(diva_fixed.elapsed / diva.elapsed, sigdigits = 3))x",
    fontsize = 12)
save("ais-momentum-fixed-vs-autotune.png", fig3)
fig3

#=
## 4. CPU vs. GPU

GPU-backed [`Chmy.Field`](@ref)s cannot be filled the way [`run_solve`](@ref) fills CPU
ones: `fill_from_grid!`/`fill_from_grid3d!` write one element at a time
(`f[i, j, k] = ...`), and CUDA.jl deliberately disallows scalar indexing on a `CuArray` — it
is almost always a performance bug, and would be one here too (each write becomes its own
kernel launch/sync). The fix is not to avoid it element-by-element; it is to never do it on
the GPU side at all: build the whole [`MechanicState`](@ref) on the CPU exactly as every
figure above already does, then move it to the GPU **in one bulk operation** —
`Adapt.adapt(CuArray, mech)`. `Field` carries its own `Adapt` rule for exactly this
(`test/api/state.jl`'s own `adapt` testset: *"a state of Fields stays adaptable — this is
what a GPU launch relies on"*), so the whole struct, every field, at its correct
location/shape, moves in one call — not a loop of per-field copies, and definitely not a
loop of per-element ones.

The mask needs the same treatment before it can be rebuilt: [`IceMask`](@ref) wraps a
`Field` reference, not a copy, so it has to be reconstructed from the *adapted* topography's
mask field, not adapted itself.

A **discarded warm-up solve** runs first and is excluded from the timed comparison:
GPU kernels compile on first use (`GPUCompiler`), and that compilation is a one-time cost of
plausibly comparable size to the ~20 s solves being measured — timing it in would answer "how
long until this GPU is warm", not "how fast does it run once warm", which is the actually
interesting number here. §2's Float64-vs-Float32 timing did **not** get this treatment (that
confound was flagged instead of fixed) because CPU JIT compilation is a much smaller fraction
of a ~20 s CPU solve than GPU kernel compilation is of a GPU one.
=#

using CUDA

if CUDA.functional()
    arch_gpu = Arch(CUDABackend())
    grid_gpu = StaggeredGrid(arch_gpu, Float64, lx, ly, dx, dy, layering)
    rt_gpu   = Runtime(grid_gpu)

    mech_cpu_for_gpu = MechanicState(grid)
    fill_from_grid!(mech_cpu_for_gpu.topography.thickness, H_ice)
    fill_from_grid!(mech_cpu_for_gpu.topography.surface, z_srf)
    fill_from_grid!(mech_cpu_for_gpu.material.viscosity_depthaveraged, visc_bar)
    fill_from_grid!(mech_cpu_for_gpu.friction.beta_eff, beta_eff)
    fill_from_grid!(mech_cpu_for_gpu.friction.beta, beta)
    fill_from_grid3d!(mech_cpu_for_gpu.material.viscosity, visc3d)
    setdata!(mech_cpu_for_gpu.velocity.depthaverage_x, 0.0)
    setdata!(mech_cpu_for_gpu.velocity.depthaverage_y, 0.0)

    topo_gpu = Pagos.Adapt.adapt(CuArray, topo)
    mask_gpu = IceMask(topo_gpu.mask.is_momentum_solved)

    function solve_on!(mech, grid, rt, mask; solver_kwargs...)
        solver = PseudoTransientSolver(grid; solver_kwargs...)
        diva_update!(mech, solver, rt, mask)
        elapsed = @elapsed result = pseudo_transient!(mech, cst, solver, rt,
                                                       DIVAMomentumBalance(), mask)
        speed = Array(Float32.(spy .* sqrt.(
            (@views (interior(mech.velocity.depthaverage_x)[1:(end - 1), :, 1] .+
                     interior(mech.velocity.depthaverage_x)[2:end, :, 1]) ./ 2) .^ 2 .+
            (@views (interior(mech.velocity.depthaverage_y)[:, 1:(end - 1), 1] .+
                     interior(mech.velocity.depthaverage_y)[:, 2:end, 1]) ./ 2) .^ 2)))
        return (; speed, elapsed, iterations = result.iterations, converged = result.converged,
               residual = result.residual)
    end

    # Warm-up: a cheap, throwaway solve (loose `abstol`, `maxiter` capped low) whose only
    # purpose is to force every kernel this solve path touches through GPUCompiler once,
    # off the clock.
    mech_warmup = Pagos.Adapt.adapt(CuArray, mech_cpu_for_gpu)
    solve_on!(mech_warmup, grid_gpu, rt_gpu, mask_gpu; SOLVER_KWARGS..., maxiter = 5, abstol = 0.0)
    mech_warmup = nothing
    GC.gc()
    println("GPU warm-up done (kernels compiled, not timed)."); ram()

    mech_gpu = Pagos.Adapt.adapt(CuArray, mech_cpu_for_gpu)
    mech_cpu_for_gpu = nothing
    GC.gc()
    diva_gpu = solve_on!(mech_gpu, grid_gpu, rt_gpu, mask_gpu; SOLVER_KWARGS...)
    mech_gpu = nothing
    GC.gc()
    println("DIVA, GPU: ", (; diva_gpu.converged, diva_gpu.iterations, diva_gpu.elapsed, diva_gpu.residual))
    ram()

    diff_gpu = on_ice(diva.speed .- diva_gpu.speed)
    diff_gpu_vals = filter(!isnan, diff_gpu)
    @printf("CPU vs GPU DIVA (on-ice): RMSE = %.4g m/yr, max|Δ| = %.4g m/yr, elapsed %.3gs vs %.3gs (%.2gx)\n",
           sqrt(mean(abs2, diff_gpu_vals)), maximum(abs, diff_gpu_vals), diva.elapsed, diva_gpu.elapsed,
           diva.elapsed / diva_gpu.elapsed)

    fig4 = Figure(size = (1150, 620))
    crange4 = (0, quantile(filter(!isnan, on_ice(diva.speed)), 0.995))
    for (col, (data, title)) in enumerate(((on_ice(diva.speed), "DIVA, CPU"),
                                           (on_ice(diva_gpu.speed), "DIVA, GPU")))
        ax = Axis(fig4[1, col], xlabel = "x (km)", ylabel = col == 1 ? "y (km)" : "",
            aspect = DataAspect(), title = title)
        hm = heatmap!(ax, xc, yc, data; colorrange = crange4)
        col == 2 && Colorbar(fig4[1, 3], hm, label = "speed (m/yr)")
    end
    Label(fig4[2, 1:3],
        "CPU: $(diva.iterations) iter, $(round(diva.elapsed, digits = 2))s   |   " *
        "GPU: $(diva_gpu.iterations) iter, $(round(diva_gpu.elapsed, digits = 2))s   |   " *
        "speedup = $(round(diva.elapsed / diva_gpu.elapsed, sigdigits = 3))x   |   " *
        "RMSE = $(round(sqrt(mean(abs2, diff_gpu_vals)), sigdigits = 3)) m/yr",
        fontsize = 12)
    save("ais-momentum-cpu-vs-gpu.png", fig4)
    fig4
else
    println("CUDA.functional() == false — skipping the CPU vs GPU comparison.")
end

#=
## 2. Pagos DIVA vs. Yelmo's own DIVA solve — surface speed, RMSE only

`ux`/`uy` are read here, not in the loading section at the top, as a **2D surface slice
straight from the file** (`ds["ux"][:, :, end, 1]`) — `zeta[end] == 1.0` is the surface (the
`zeta` printed while writing this script runs bed-to-surface,
`[0.0, 0.01, ..., 1.0]`). The other 10 layers of the full 3D `ux`/`uy` are never read into
memory at all, which is the whole point: this comparison, per its own scoping decision, is
surface-only, so there is nothing to gain from the other layers and no reason to pay for
them. Float16 is safe here — Yelmo's surface speeds top out under 4300 m/yr in this file,
comfortably inside Float16's ~65504 magnitude cap (checked, like every other narrow-type
load above, against the file's own values rather than assumed).
=#

ux_srf, uy_srf = NCDataset(restart_file) do ds
    (Float16.(ds["ux"][:, :, end, 1]), Float16.(ds["uy"][:, :, end, 1]))
end
speed_yelmo_srf = on_ice(sqrt.(Float32.(ux_srf) .^ 2 .+ Float32.(uy_srf) .^ 2))

# `diva.speed` (§1) is the *depth-averaged* speed, not the surface one — comparing it
# against Yelmo's surface speed would be comparing two different physical quantities.
# `velocities3D!`'s Eq. 17 surface reconstruction is what belongs on this side of the
# comparison; §1 only ran the depth-averaged solve, not the reconstruction, so this section
# reruns DIVA (reusing the same warm kernels from §1) and reconstructs the surface velocity
# from the converged state before comparing.
mech_srf = MechanicState(grid)
fill_from_grid!(mech_srf.topography.thickness, H_ice)
fill_from_grid!(mech_srf.topography.surface, z_srf)
fill_from_grid!(mech_srf.material.viscosity_depthaveraged, visc_bar)
fill_from_grid!(mech_srf.friction.beta_eff, beta_eff)
fill_from_grid!(mech_srf.friction.beta, beta)
fill_from_grid3d!(mech_srf.material.viscosity, visc3d)
setdata!(mech_srf.velocity.depthaverage_x, 0.0)
setdata!(mech_srf.velocity.depthaverage_y, 0.0)

solver_srf = PseudoTransientSolver(grid; SOLVER_KWARGS...)
diva_update!(mech_srf, solver_srf, rt, mask)
pseudo_transient!(mech_srf, cst, solver_srf, rt, DIVAMomentumBalance(), mask)
velocities3D!(mech_srf, rt, DIVAMomentumBalance(), mask)

# Staggered onto cell centres by averaging adjacent faces — not a truncation — matching
# the exact pattern `run_solve` already uses for `depthaverage_x`/`y`. `surface_x` is
# `ACX2` (`nx+1, ny`), `surface_y` is `ACY2` (`nx, ny+1`); averaging each along its own
# Vertex axis gives both the same `(nx, ny)` cell-centred shape.
speed_pagos_srf = on_ice(Float32.(spy .* sqrt.(
    (@views (interior(mech_srf.velocity.surface_x)[1:(end - 1), :, 1] .+
             interior(mech_srf.velocity.surface_x)[2:end, :, 1]) ./ 2) .^ 2 .+
    (@views (interior(mech_srf.velocity.surface_y)[:, 1:(end - 1), 1] .+
             interior(mech_srf.velocity.surface_y)[:, 2:end, 1]) ./ 2) .^ 2)))
mech_srf = nothing
GC.gc()

diff_srf_vals = filter(!isnan, speed_pagos_srf .- speed_yelmo_srf)
@printf("Pagos DIVA vs Yelmo, surface speed (on-ice, %d cells): RMSE = %.4g m/yr, mean|Δ| = %.4g m/yr, max|Δ| = %.4g m/yr\n",
       length(diff_srf_vals), sqrt(mean(abs2, diff_srf_vals)), mean(abs, diff_srf_vals),
       maximum(abs, diff_srf_vals))
ram()

fig5 = Figure(size = (1150, 620))
crange5 = (0, quantile(filter(!isnan, speed_yelmo_srf), 0.995))
for (col, (data, title)) in enumerate(((speed_pagos_srf, "Pagos DIVA (surface)"),
                                       (speed_yelmo_srf, "Yelmo (surface)")))
    ax = Axis(fig5[1, col], xlabel = "x (km)", ylabel = col == 1 ? "y (km)" : "",
        aspect = DataAspect(), title = title)
    hm = heatmap!(ax, xc, yc, data; colorrange = crange5)
    col == 2 && Colorbar(fig5[1, 3], hm, label = "speed (m/yr)")
end
Label(fig5[2, 1:3],
    "RMSE = $(round(sqrt(mean(abs2, diff_srf_vals)), sigdigits = 3)) m/yr   |   " *
    "mean|Δ| = $(round(mean(abs, diff_srf_vals), sigdigits = 3)) m/yr   |   " *
    "max|Δ| = $(round(maximum(abs, diff_srf_vals), sigdigits = 3)) m/yr   |   " *
    "n = $(length(diff_srf_vals)) cells",
    fontsize = 12)
save("ais-momentum-pagos-vs-yelmo.png", fig5)
fig5

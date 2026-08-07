#=

# Performance scaling: linear solve & Gershgorin-tuned PT vs. resolution

Wall-clock sweep of four DIVA momentum-solve configurations over 32/16/8/4 km:

 1. **Linear, CPU** — direct sparse solve ([`LinearMomentumSolver2D`](@ref), `SparseArrays.lu`).
 2. **Linear, GPU** — the same solver's CUDSS-backed `velocity!` (`ext/PagosCUDSSExt.jl`).
 3. **Gershgorin PT, CPU** — [`PseudoTransientSolver`](@ref), library defaults.
 4. **Gershgorin PT, GPU** — same, on `CUDABackend()`.

Timings land in `perf-scaling.csv` (via `DelimitedFiles`); `perf-scaling-lines.jl` reads that
file and plots wall time against resolution. One sanity figure per resolution is saved
alongside, so a leg that "finished fast" because it solved the wrong problem is visible
rather than buried in a number.

## Why every resolution comes from the 8 km restart

A real 16 km Yelmo restart exists on disk (see `helpers.jl`), but using it would confound two
things this sweep needs kept apart: how solve time scales with resolution, and how Yelmo's own
state happens to differ between its restart runs. Every resolution here is instead derived
from the *same* 8 km restart, so grid spacing is the only thing that changes between legs.

`coarsen` is a non-overlapping block mean, used for 32 and 16 km; 8 km is the restart passed
through untouched; `refine` is bilinear interpolation, used for 4 km. Neither is a
conservative remap:

  - Coarsening turns `f_grnd` into a sub-grid grounded *fraction* (harmless — `is_grounded`
    thresholds it at `> 0` anyway) and can dilute thin or isolated ice below the `H_min = 0`
    `is_ice` threshold. At 32 km this already erases all 48 detached iceberg cells the 8 km
    restart carries, so `n_detached` drops to 0 — expected, not a bug.
  - Refining adds no information. The 4 km leg resolves nothing the 8 km field did not
    already contain, and bilinear interpolation smears `H` across the ice margin. It is a
    **size stress-test — a genuinely 4× larger, well-posed system to time — not a 4 km
    simulation**, and its wall times should be read that way.

## Convergence tolerance: why `abstol = 2e-3` here and not `1e-3`

[`ScaledResidual`](@ref) measures a *rate*: `_driving_rate!` forms `τ_d/(ρH)` and `err` is
its max-norm over the domain-max scale. The `1/H` factor means **the thinnest ice in the
domain sets the norm**, so `err` is decided by a handful of margin cells rather than by the
bulk of the ice sheet — and which cells those are moves around with resolution.

Measured on this sweep at `abstol = 1e-3`: 16 km stalls at `err ≈ 1.43e-3` from exactly
**3 faces out of 144 780** (0.002%), one patch of thin (H = 27–66 m), floating
(`β_eff = 0`) margin ice, and burns all 2000 iterations — while its velocity field is
indistinguishable from the direct linear solve. 8 km clears the same bar at `0.81e-3`, i.e.
by only 20%. The result is a non-monotonic iteration count (140 → 2000 → 340) that reflects
where the thinnest margin cell happens to land, not solver cost.

Three things that do **not** fix it, all tested: raising `icemasks!`'s `H_min` (0/1/5/10 m —
no effect, the offending cells are 27–66 m thick); aggregating with a minimum ice fraction
per coarse cell (0.25/0.5 — 3–7× *worse*, since dropping cells sharpens margins and isolates
more thin patches); and switching to [`VelocityIncrement`](@ref), which the
[`pseudo_transient!`](@ref) docstring rules out for exactly this kind of domain ("the
increment criterion cannot distinguish a converged shelf from a stalled one").

So `abstol = 2e-3` sits just above that thin-margin floor, which is what makes the four legs
comparable: every one of them then reports genuine time-to-solution. It is a deliberate
departure from the `1e-3` the other scripts here use, and it makes PT look slightly better
against the direct solve than a tighter tolerance would.

## Timing convention

`diva_update!`'s β_eff correction runs once per resolution *outside* every timed region,
matching `run_solve` (`helpers.jl`) and the linear-solve scripts. All four legs therefore time
only the core solve on an identical, already-DIVA-corrected problem: `A_active \\ b_active` /
CUDSS for the linear legs, `pseudo_transient!` for the PT legs. Each GPU leg runs a discarded
warm-up first (CUDSS symbolic analysis; PT kernel compilation past `SOLVER_KWARGS`' tuning
cadence), and the **whole 32 km resolution is run twice with the first pass discarded**, so
CPU-side JIT of the assembly loops and both solver paths is paid before anything is recorded.

!!! warning "Memory: 4 km needs a mostly-idle machine"
    The 4 km leg holds ~3.9 GiB of host **and** ~3.9 GiB of device memory for its
    [`MechanicState`](@ref) simultaneously (`Adapt.adapt` copies rather than moves), on top of
    ~0.7 GiB of sparse operator and the `lu` factors. Every leg drops its state and forces a
    `GC.gc()` (plus `CUDA.reclaim()`) before the next begins, so peak usage is one resolution
    at a time — but 4 km still wants ~8 GiB of free RAM and a near-empty 8 GiB GPU.

    Watch **swap** specifically, not just free RAM: with swap exhausted, the 8 km `lu`
    degraded into thrashing that ran >100× slower rather than failing outright, which is much
    harder to recognise than an `OutOfMemoryError`.
=#
using Pkg
Pkg.activate(joinpath(@__DIR__, "../../.."))
using CUDA
if CUDA.functional()
    CUDA.zeros(1)
    CUDA.synchronize()
end
using Pagos, NCDatasets, CairoMakie, Statistics, Printf
using SparseArrays, LinearAlgebra, CUDA.CUSPARSE, CUDSS, DelimitedFiles

figdir = joinpath(@__DIR__, "figs")
set_theme!(theme_latexfonts())

#=
## Carrying speed fields between split processes

A resolution split across processes (4 km) never holds all four solutions at once, so each
pass writes the quadrants it computed as raw `Float32` and the figure is rendered as soon as
every quadrant is on disk. Raw dumps rather than a serialization format: these are plain
dense `nx × ny` arrays whose shape is already known from the grid, so there is nothing to
encode.
=#
speeddir = joinpath(figdir, ".speeds")
speedpath(km, name) = joinpath(speeddir, "$(km)km-$(name).bin")

# Each leg persists its speed field plus the three scalars its panel subtitle needs, so a
# figure assembled in one process can label a leg that ran in the other.
function dump_leg(km, name, speed, elapsed, iters, converged)
    mkpath(speeddir)
    write(speedpath(km, name), Float32.(speed))
    open(speedpath(km, name) * ".meta", "w") do io
        println(io, elapsed); println(io, iters); println(io, converged)
    end
    return nothing
end

function load_leg(km, name, nx, ny)
    p = speedpath(km, name)
    isfile(p) || return (fill(NaN32, nx, ny), NaN, 0, false)
    speed = reshape(reinterpret(Float32, read(p)), nx, ny)
    m = isfile(p * ".meta") ? readlines(p * ".meta") : ["NaN", "0", "false"]
    return (speed, parse(Float64, m[1]), parse(Int, m[2]), parse(Bool, m[3]))
end

hasdata(a) = any(isfinite, a)

#=
## Load the 8 km restart once — every resolution below is resampled from these fields
=#
restart_file = "/home/jan/pCloudSync/PhD/Projects/Ice-Sheet-Modelling/ice-data-pagos/yelmo_restart_ais_8km.nc"
load2d(ds, name, T) = T.(dropdims(ds[name][:, :, :]; dims = 3))
load3d(ds, name, T) = T.(dropdims(ds[name][:, :, :, :]; dims = 4))

xc0, yc0, H_ice0, z_srf0, z_bed0, f_grnd0, visc_bar0, beta0, beta_eff0, visc3d0, zeta =
NCDataset(restart_file) do ds
    (Float64.(ds["xc"][:]), Float64.(ds["yc"][:]),
     load2d(ds, "H_ice", Float32), load2d(ds, "z_srf", Float32), load2d(ds, "z_bed", Float32),
     load2d(ds, "f_grnd", Float32),
     load2d(ds, "visc_bar", Float32), load2d(ds, "beta", Float32), load2d(ds, "beta_eff", Float32),
     load3d(ds, "visc", Float32),
     Float64.(ds["zeta"][:]))
end
nz  = length(zeta)
dx0 = (xc0[2] - xc0[1]) * 1e3          # m
dxkm0 = xc0[2] - xc0[1]                # km, for plot axes

#=
## Resampling

`ratio = target_dx / 8 km`: `> 1` coarsens by block mean, `< 1` refines by bilinear
interpolation, `== 1` passes the field through untouched.
=#
function coarsen(data::AbstractMatrix, factor::Int)
    ni, nj = size(data)
    no, mo = ni ÷ factor, nj ÷ factor
    out = zeros(eltype(data), no, mo)
    @inbounds for jo in 1:mo, io in 1:no
        acc = zero(eltype(data))
        for dj in 1:factor, di in 1:factor
            acc += data[(io - 1) * factor + di, (jo - 1) * factor + dj]
        end
        out[io, jo] = acc / factor^2
    end
    return out
end

# Bilinear onto the refined cell centres: refined cell `io` sits at fractional old-index
# `(io - 0.5)/factor + 0.5`, i.e. the two children of old cell `i` land at `i ∓ 0.25`.
function refine(data::AbstractMatrix, factor::Int)
    ni, nj = size(data)
    no, mo = ni * factor, nj * factor
    out = zeros(eltype(data), no, mo)
    @inbounds for jo in 1:mo, io in 1:no
        xi = (io - 0.5) / factor + 0.5
        yj = (jo - 0.5) / factor + 0.5
        i0 = clamp(floor(Int, xi), 1, ni - 1)
        j0 = clamp(floor(Int, yj), 1, nj - 1)
        tx = clamp(xi - i0, 0, 1)
        ty = clamp(yj - j0, 0, 1)
        out[io, jo] = (1 - tx) * (1 - ty) * data[i0, j0] + tx * (1 - ty) * data[i0 + 1, j0] +
                      (1 - tx) * ty * data[i0, j0 + 1] + tx * ty * data[i0 + 1, j0 + 1]
    end
    return out
end

function resample(data::AbstractMatrix, ratio::Rational)
    ratio == 1 && return copy(data)
    ratio > 1 && return coarsen(data, Int(ratio))
    return refine(data, Int(inv(ratio)))
end
resample(data::AbstractArray{<:Any,3}, ratio::Rational) =
    cat((resample(view(data, :, :, k), ratio) for k in axes(data, 3))...; dims = 3)

#=
## Per-resolution helpers

Copied from `helpers.jl` rather than `include`d: that file hardcodes one resolution's restart
file and builds a single grid at module scope, which is exactly what this sweep must not do.
=#
function fill_from_grid!(f, data)
    ni, nj = size(data)
    for k in axes(interior(f), 3), j in -1:(nj + 2), i in -1:(ni + 2)
        ic, jc = clamp(i, 1, ni), clamp(j, 1, nj)
        f[i, j, k] = data[ic, jc]
    end
    return f
end

function fill_from_grid3d!(f, data)
    ni, nj, nk = size(data)
    for k in 1:nk, j in -1:(nj + 2), i in -1:(ni + 2)
        ic, jc = clamp(i, 1, ni), clamp(j, 1, nj)
        f[i, j, k] = data[ic, jc, k]
    end
    return f
end

#=
`abstol = 2e-3`, not the `1e-3` every other script in this folder uses — see the
"Convergence tolerance" note in the module docstring for why 1e-3 is inside the noise band
of this norm on real margin geometry.
=#
const SOLVER_KWARGS = (
    abstol = 2e-3, maxiter = 2000, ncheck = 20, printout_every = 50,
    pseudo_timestep = GershgorinPseudoTimeStep(cfl = 0.8),
    convergence = ScaledResidual(),
    friction_update = ActiveFrictionUpdate(),
    tuning = AutotunedDynamicRelaxation(cadence = 50),
)
const CST = Constants{Float32}()

#=
Leg selection, so a resolution whose state does not fit alongside everything else can be
split across processes. At 4 km `build_mech` alone takes RSS to ~9.4 GiB, which leaves no
room to also hold the device copy, so its CPU and GPU legs must run separately:

```
PAGOS_SKIP_GPU=1 julia --project=docs .../perf-scaling.jl 4   # PT CPU + linear CPU
PAGOS_ONLY_GPU=1 julia --project=docs .../perf-scaling.jl 4   # PT GPU + linear GPU
```

A leg that does not run records `NaN`; the CSV writer merges those fragments into the single
row for that resolution, so the plotting script still sees one row per resolution.
=#
const RUN_CPU_LEGS = get(ENV, "PAGOS_ONLY_GPU", "0") != "1"
const RUN_GPU_LEGS = get(ENV, "PAGOS_SKIP_GPU", "0") != "1"

depthavg_speed(mech) = Float32.(sqrt.(
    (@views (interior(mech.velocity.depthaverage_x)[1:(end - 1), :, 1] .+
             interior(mech.velocity.depthaverage_x)[2:end, :, 1]) ./ 2) .^ 2 .+
    (@views (interior(mech.velocity.depthaverage_y)[:, 1:(end - 1), 1] .+
             interior(mech.velocity.depthaverage_y)[:, 2:end, 1]) ./ 2) .^ 2))

#=
## One resolution: build the state, run all four solves, time each core solve
=#
# Build the DIVA-corrected mechanic state. Kept out of `run_resolution` so it does not
# capture that function's locals in a closure — at 4 km a stray capture of a ~6.5 GiB state
# is the difference between fitting in RAM and being OOM-killed.
function build_mech(grid, rt, mask, H_ice, z_srf, visc_bar, beta, beta_eff, visc3d, ::Type{T}) where {T}
    m = MechanicState(grid)
    fill_from_grid!(m.topography.thickness, H_ice)
    fill_from_grid!(m.topography.surface, z_srf)
    fill_from_grid!(m.material.viscosity_depthaveraged, visc_bar)
    fill_from_grid!(m.friction.beta_eff, beta_eff)
    fill_from_grid!(m.friction.beta, beta)
    fill_from_grid3d!(m.material.viscosity, visc3d)
    setdata!(m.velocity.depthaverage_x, zero(T))
    setdata!(m.velocity.depthaverage_y, zero(T))
    s = PseudoTransientSolver(grid; SOLVER_KWARGS...)
    diva_update!(m, s, rt, mask)   # β_eff correction — excluded from every timing
    return m, s
end

function run_resolution(ratio::Rational; save_figure = true)
    T = Float32
    dx = dy = dx0 * ratio
    resolution_km = round(Int, dx / 1000)
    @printf("\n=== %d km (ratio = %s) ===\n", resolution_km, ratio)

    H_ice    = resample(H_ice0, ratio)
    z_srf    = resample(z_srf0, ratio)
    z_bed    = resample(z_bed0, ratio)
    visc_bar = resample(visc_bar0, ratio)
    beta     = resample(beta0, ratio)
    beta_eff = resample(beta_eff0, ratio)
    visc3d   = resample(visc3d0, ratio)
    f_grnd   = resample(f_grnd0, ratio)
    nx, ny   = size(H_ice)
    lx, ly   = nx * dx, ny * dy

    visc_off_ice   = maximum(visc_bar)
    visc3d_off_ice = maximum(visc3d)
    @. visc_bar = ifelse(H_ice > 0, visc_bar, visc_off_ice)
    @. visc3d   = ifelse(H_ice > 0, visc3d, visc3d_off_ice)
    on_ice(a) = ifelse.(H_ice .> 0, a, NaN)

    layering = CorrectedVerticalLayering(T, QuadraticSigmaTransform(T, nz))
    grid  = StaggeredGrid(T, lx, ly, dx, dy, layering)
    rt    = Runtime(grid)
    masks = TopographyMasks(grid)

    thickness_ice = Field(grid.arch, grid.grid2d, (Center(), Center(), Center()), T; halo = 1)
    fill_from_grid!(thickness_ice, H_ice)
    fill_from_grid!(masks.is_grounded, f_grnd .> 0)
    icemasks!(masks, thickness_ice, rt)
    momentum_mask!(masks.is_momentum_solved, masks.is_ice, masks.is_grounded, rt)
    mask = IceMask(masks.is_momentum_solved)
    n_detached = count(asarray(masks.is_ice) .& .!asarray(masks.is_momentum_solved))

    #=
    ### Peak-memory staging, and why PT runs before the linear legs

    At 4 km the [`MechanicState`](@ref) is ~6.5 GiB live (the 38 column fields carry a halo in
    `z` as well as in `x`/`y`), the assembled `LinearMomentumSolver2D` ~2 GiB, and the `lu`
    factors ~2.5 GiB — but none of them need to coexist, so the order below is chosen to keep
    them apart rather than to read naturally:

      1. `mech` is built **once**. An earlier version built it a second time for the PT legs
         after the factorization; measured, that was the actual peak (RSS 9.4 GiB, killed by
         `earlyoom`), because Julia does not return the first copy's pages to the OS.
      2. The four fields the linear operator needs are snapshotted out of `mech` (~50 MB)
         *before* PT touches the velocity, so the linear legs can run long after `mech` is gone.
      3. PT runs first, and the host state is dropped the moment it has been copied to the
         device — so the GPU leg holds ~6.5 GiB on the card and almost nothing on the host.
      4. The linear legs run last, on the snapshot, with `mech` already released.
    =#
    mech, solver_cpu = build_mech(grid, rt, mask, H_ice, z_srf, visc_bar, beta, beta_eff, visc3d, T)

    # Snapshot for the linear operator (step 2 above): copies, not views into `mech`.
    Hf    = copy(interior(mech.topography.thickness)[:, :, 1])
    sf    = copy(interior(mech.topography.surface)[:, :, 1])
    viscf = copy(interior(mech.material.viscosity_depthaveraged)[:, :, 1])
    befff = copy(interior(mech.friction.beta_eff)[:, :, 1])
    run_gpu = CUDA.functional() && RUN_GPU_LEGS

    #=
    ### 1. Gershgorin PT, CPU
    =#
    if RUN_CPU_LEGS
        elapsed_pt_cpu = @elapsed res_pt_cpu =
            pseudo_transient!(mech, CST, solver_cpu, rt, DIVAMomentumBalance(), mask)
        speed_pt_cpu = on_ice(depthavg_speed(mech))
        @printf("  PT CPU:     %8.2f s  (%d it, converged = %s)\n",
               elapsed_pt_cpu, res_pt_cpu.iterations, res_pt_cpu.converged)
    else
        elapsed_pt_cpu = NaN
        res_pt_cpu = (; converged = false, iterations = 0)
        speed_pt_cpu = fill(NaN32, nx, ny)
    end

    #=
    ### 2. Gershgorin PT, GPU
    =#
    if run_gpu
        setdata!(mech.velocity.depthaverage_x, zero(T))
        setdata!(mech.velocity.depthaverage_y, zero(T))

        arch_gpu  = Arch(CUDABackend())
        grid_gpu  = StaggeredGrid(arch_gpu, T, lx, ly, dx, dy, layering)
        rt_gpu    = Runtime(grid_gpu)
        masks_gpu = Pagos.Adapt.adapt(CuArray, masks)
        mask_gpu  = IceMask(masks_gpu.is_momentum_solved)
        mech_gpu  = Pagos.Adapt.adapt(CuArray, mech)
        # Host state is dead once it is on the card: releasing it here is what keeps the
        # 4 km GPU leg from holding ~6.5 GiB on both sides at the same time.
        mech = nothing; solver_cpu = nothing
        GC.gc()

        function solve_on!(m, g, r, msk; solver_kwargs...)
            solver = PseudoTransientSolver(g; solver_kwargs...)
            elapsed = @elapsed result =
                pseudo_transient!(m, CST, solver, r, DIVAMomentumBalance(), msk)
            return elapsed, result
        end

        # Warm-up past SOLVER_KWARGS' tuning cadence (50) so `_tune!` compiles off the clock.
        solve_on!(mech_gpu, grid_gpu, rt_gpu, mask_gpu;
            (; SOLVER_KWARGS..., maxiter = 55, abstol = 0.0)...)
        setdata!(mech_gpu.velocity.depthaverage_x, zero(T))
        setdata!(mech_gpu.velocity.depthaverage_y, zero(T))

        elapsed_pt_gpu, res_pt_gpu = solve_on!(mech_gpu, grid_gpu, rt_gpu, mask_gpu; SOLVER_KWARGS...)
        speed_pt_gpu = on_ice(Array(depthavg_speed(mech_gpu)))
        @printf("  PT GPU:     %8.2f s  (%d it, converged = %s)\n",
               elapsed_pt_gpu, res_pt_gpu.iterations, res_pt_gpu.converged)

        mech_gpu = nothing; masks_gpu = nothing; mask_gpu = nothing
        GC.gc(); CUDA.reclaim()
    else
        mech = nothing; solver_cpu = nothing
        GC.gc()
        elapsed_pt_gpu = NaN
        res_pt_gpu = (; converged = false, iterations = 0)
        speed_pt_gpu = fill(NaN32, nx, ny)
    end

    i_idx = PeriodicIndexing(1, nx)
    j_idx = PeriodicIndexing(1, ny)
    rho_ice, grav = T(CST.density_ice), T(CST.gravity)

    N     = similar(Hf); N_ab = similar(Hf)
    b_acx = similar(Hf); b_acy = similar(Hf)
    tx    = similar(Hf); ty   = similar(Hf)
    for i in 1:nx, j in 1:ny
        im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
        N[i, j] = Hf[i, j] * viscf[i, j]
        eta_ab = 4 / (inv(viscf[i, j]) + inv(viscf[ip1, j]) + inv(viscf[i, jp1]) + inv(viscf[ip1, jp1]))
        H_ab   = (Hf[i, j] + Hf[ip1, j] + Hf[i, jp1] + Hf[ip1, jp1]) / 4
        N_ab[i, j] = eta_ab * H_ab
        b_acx[i, j] = (befff[i, j] + befff[ip1, j]) / 2
        b_acy[i, j] = (befff[i, j] + befff[i, jp1]) / 2
        Hx = (Hf[i, j] + Hf[ip1, j]) / 2
        Hy = (Hf[i, j] + Hf[i, jp1]) / 2
        tx[i, j] = rho_ice * grav * Hx * (sf[ip1, j] - sf[i, j]) / dx
        ty[i, j] = rho_ice * grav * Hy * (sf[i, jp1] - sf[i, j]) / dy
    end

    regular_grid = RegularGrid(T, (nx - 1) * dx, (ny - 1) * dy, dx, dy)
    lsd = LinearMomentumSolver2D(regular_grid, DIVAMomentumBalance())
    ux0 = zeros(T, nx, ny); uy0 = zeros(T, nx, ny)
    populate_vectors!(lsd, N, N_ab, ux0, uy0, tx, ty, b_acx, b_acy)

    is_solved = interior(masks.is_momentum_solved)[:, :, 1]
    active_dof = falses(2 * nx * ny)
    for i in 1:nx, j in 1:ny
        im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
        active_dof[Pagos._ij2n_ux(i, j, nx, ny)] = is_solved[i, j] | is_solved[ip1, j]
        active_dof[Pagos._ij2n_uy(i, j, nx, ny)] = is_solved[i, j] | is_solved[i, jp1]
    end
    idx = findall(active_dof)
    A_active = lsd.A[idx, idx]
    b_active = lsd.b[idx]
    lsd = nothing          # ~2 GiB of pattern/permutation arrays at 4 km, dead from here
    GC.gc()

    # `velocity!(ux, uy, lsd)`'s index mapping, inlined so `lsd` can be dropped above rather
    # than kept alive across the factorization just to unpack the solution vector.
    function extract_speed(u_solved)
        u = zeros(T, 2 * nx * ny)
        u[idx] .= T.(u_solved)
        ux = zeros(T, nx, ny); uy = zeros(T, nx, ny)
        for i in 1:nx, j in 1:ny
            ux[i, j] = u[Pagos._ij2n_ux(i, j, nx, ny)]
            uy[i, j] = u[Pagos._ij2n_uy(i, j, nx, ny)]
        end
        ux_c = similar(ux); uy_c = similar(uy)
        for i in 1:nx, j in 1:ny
            im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
            ux_c[i, j] = (ux[im1, j] + ux[i, j]) / 2
            uy_c[i, j] = (uy[i, jm1] + uy[i, j]) / 2
        end
        return T.(sqrt.(ux_c .^ 2 .+ uy_c .^ 2))
    end

    #=
    ### 3. Linear, CPU
    =#
    if RUN_CPU_LEGS
        elapsed_lin_cpu = @elapsed u_cpu = A_active \ b_active
        speed_lin_cpu = on_ice(extract_speed(u_cpu))
        u_cpu = nothing; GC.gc()
        @printf("  linear CPU: %8.2f s  (%d active DOFs)\n", elapsed_lin_cpu, length(idx))
    else
        elapsed_lin_cpu = NaN
        speed_lin_cpu = fill(NaN32, nx, ny)
    end

    #=
    ### 4. Linear, GPU (CUDSS)

    `lsd_gpu` is dropped and the device pool reclaimed after each call, so the CUDSS factors
    do not outlive the solve that needed them. By this point the PT legs have already run and
    released `mech_gpu`, so the card is otherwise empty.
    =#
    function gpu_cudss_solve(A_active, b_active)
        A_gpu = CuSparseMatrixCSR(A_active)
        b_gpu = CuArray(b_active)
        u_gpu = CUDA.zeros(T, length(b_active))
        lsd_gpu = LinearMomentumSolver2D(
            DIVAMomentumBalance(), 0, 0, zero(T), zero(T), zero(T),
            u_gpu, similar(u_gpu), b_gpu, A_gpu, CUDA.zeros(Int, 1),
            PeriodicIndexing(1, 1), PeriodicIndexing(1, 1), Ref{Any}(nothing),
        )
        velocity!(lsd_gpu)
        CUDA.synchronize()
        return Array(lsd_gpu.u)
    end

    if run_gpu
        gpu_cudss_solve(A_active, b_active)              # warm-up (CUDSS symbolic analysis)
        GC.gc(); CUDA.reclaim()
        elapsed_lin_gpu = @elapsed u_gpu = gpu_cudss_solve(A_active, b_active)
        speed_lin_gpu = on_ice(extract_speed(u_gpu))
        GC.gc(); CUDA.reclaim()
        @printf("  linear GPU: %8.2f s\n", elapsed_lin_gpu)
    else
        elapsed_lin_gpu = NaN
        speed_lin_gpu = fill(NaN32, nx, ny)
    end

    A_active = nothing; b_active = nothing
    GC.gc()

    #=
    Persist whatever this pass solved, then reload all four quadrants — a leg run in a sibling
    process comes back off disk, one skipped everywhere comes back as `NaN`. The figure is
    rendered only when all four are present, so a split run does not clobber a complete
    figure with half-empty panels.
    =#
    if RUN_CPU_LEGS
        dump_leg(resolution_km, "pt_cpu", speed_pt_cpu, elapsed_pt_cpu,
                 res_pt_cpu.iterations, res_pt_cpu.converged)
        dump_leg(resolution_km, "lin_cpu", speed_lin_cpu, elapsed_lin_cpu, 0, true)
    end
    if run_gpu
        dump_leg(resolution_km, "pt_gpu", speed_pt_gpu, elapsed_pt_gpu,
                 res_pt_gpu.iterations, res_pt_gpu.converged)
        dump_leg(resolution_km, "lin_gpu", speed_lin_gpu, elapsed_lin_gpu, 0, true)
    end

    s_pt_cpu,  e_pt_cpu,  it_pt_cpu,  cv_pt_cpu  = load_leg(resolution_km, "pt_cpu", nx, ny)
    s_pt_gpu,  e_pt_gpu,  it_pt_gpu,  cv_pt_gpu  = load_leg(resolution_km, "pt_gpu", nx, ny)
    s_lin_cpu, e_lin_cpu, _, _                   = load_leg(resolution_km, "lin_cpu", nx, ny)
    s_lin_gpu, e_lin_gpu, _, _                   = load_leg(resolution_km, "lin_gpu", nx, ny)

    if save_figure && all(hasdata, (s_pt_cpu, s_pt_gpu, s_lin_cpu, s_lin_gpu))
        save_resolution_figure(; resolution_km, ratio, nx, ny, z_bed,
            elapsed_lin_cpu = e_lin_cpu, elapsed_lin_gpu = e_lin_gpu,
            elapsed_pt_cpu = e_pt_cpu, elapsed_pt_gpu = e_pt_gpu,
            res_pt_cpu = (; iterations = it_pt_cpu, converged = cv_pt_cpu),
            res_pt_gpu = (; iterations = it_pt_gpu, converged = cv_pt_gpu),
            run_gpu = true,
            speed_lin_cpu = s_lin_cpu, speed_lin_gpu = s_lin_gpu,
            speed_pt_cpu = s_pt_cpu, speed_pt_gpu = s_pt_gpu)
    end

    out = (; resolution_km, nx, ny, n_dof = length(idx), n_detached,
          elapsed_lin_cpu, elapsed_lin_gpu, elapsed_pt_cpu, elapsed_pt_gpu,
          pt_cpu_iterations = res_pt_cpu.iterations, pt_gpu_iterations = res_pt_gpu.iterations,
          pt_cpu_converged = res_pt_cpu.converged, pt_gpu_converged = res_pt_gpu.converged)

    masks = nothing
    GC.gc()
    return out
end

#=
## Sanity figure, one per resolution

Crop scales with resolution so the same physical margin is trimmed at every grid size.
=#
function save_resolution_figure(; resolution_km, ratio, nx, ny, z_bed,
        elapsed_lin_cpu, elapsed_lin_gpu, elapsed_pt_cpu, elapsed_pt_gpu,
        res_pt_cpu, res_pt_gpu, run_gpu,
        speed_lin_cpu, speed_lin_gpu, speed_pt_cpu, speed_pt_gpu)

    crop_x = clamp(round(Int, 30 / ratio), 1, nx ÷ 2 - 1)
    crop_y = clamp(round(Int, 80 / ratio), 1, ny ÷ 2 - 1)
    ix = (crop_x + 1):(nx - crop_x)
    iy = (crop_y + 1):(ny - crop_y)
    dxkm = dxkm0 * ratio
    xc_c = [xc0[1] - dxkm0 / 2 + (i - 0.5) * dxkm for i in ix]
    yc_c = [yc0[1] - dxkm0 / 2 + (j - 0.5) * dxkm for j in iy]
    z_bed_c = z_bed[ix, iy]
    crop(d) = d[ix, iy]

    stops = [0, 20, 100, 400, 700, 1000]
    speed_cmap = cgrad(
        [:white, :white, :dodgerblue4, :lightgoldenrod1, :orangered, :darkred],
        stops ./ 1000,
    )
    crange = extrema(stops)

    function panel!(fig_pos, title, subtitle, data)
        ax = Axis(fig_pos, title = title, subtitle = subtitle, aspect = DataAspect())
        heatmap!(ax, xc_c, yc_c, z_bed_c; colormap = :oleron, colorrange = (-6000, 6000))
        hm = heatmap!(ax, xc_c, yc_c, max.(data, 1); colorrange = crange,
            colormap = speed_cmap, lowclip = speed_cmap[1], highclip = speed_cmap[end])
        hidedecorations!(ax)
        return hm
    end

    fig = Figure(size = (900, 950), fontsize = 18)
    Label(fig[0, 1:2], "$(resolution_km) km  ($(nx)×$(ny))", fontsize = 22, font = :bold)
    panel!(fig[1, 1], "Linear, CPU", @sprintf("%.2f s", elapsed_lin_cpu), crop(speed_lin_cpu))
    hm = panel!(fig[1, 2], "Linear, GPU (CUDSS)",
        run_gpu ? @sprintf("%.2f s", elapsed_lin_gpu) : "CUDA unavailable", crop(speed_lin_gpu))
    panel!(fig[2, 1], "Gershgorin PT, CPU",
        @sprintf("%d it · %.2f s", res_pt_cpu.iterations, elapsed_pt_cpu), crop(speed_pt_cpu))
    panel!(fig[2, 2], "Gershgorin PT, GPU",
        run_gpu ? @sprintf("%d it · %.2f s", res_pt_gpu.iterations, elapsed_pt_gpu) : "CUDA unavailable",
        crop(speed_pt_gpu))
    Colorbar(fig[3, 1:2], hm, vertical = false, width = Relative(0.5), flipaxis = false,
        height = Relative(1), label = "speed (m/yr)")
    rowsize!(fig.layout, 3, 20)
    colgap!(fig.layout, 5)
    save("$figdir/perf-scaling-$(resolution_km)km.png", fig)
    return nothing
end

#=
## Run the sweep

Coarse → fine, so a memory-driven failure at the finest leg still leaves every cheaper one
recorded. The 32 km leg runs twice and the first pass is discarded: it pays the CPU-side JIT
for the assembly loops, both solver paths and the plotting stack, which would otherwise land
entirely inside the first recorded resolution.

**4 km must run in its own process**, which is why it is not in the default set:

```
julia --project=docs docs/src/examples/ais-momentum/perf-scaling.jl        # 32, 16, 8 km
julia --project=docs docs/src/examples/ais-momentum/perf-scaling.jl 4      # 4 km, appends
```

Not a stylistic preference — measured. Julia does not return a freed `MechanicState`'s pages
to the OS, so after the 8 km leg the process already carries several GiB of heap high-water;
4 km's ~6.5 GiB state then pushes RSS to ~10.9 GiB and `earlyoom` kills it. The same 4 km leg
in a *fresh* process peaks at ~9.4 GiB and completes. The 4 km invocation also skips the
32 km warm-up pass, since that warm-up's own high-water is enough to reintroduce the problem;
the JIT it would have paid is a few seconds against 4 km solves of tens to hundreds, so the
distortion is far smaller than the risk.
=#
const RATIO_OF_KM = Dict(32 => 4//1, 16 => 2//1, 8 => 1//1, 4 => 1//2)

requested_km = isempty(ARGS) ? [32, 16, 8] : parse.(Int, ARGS)
default_run  = isempty(ARGS)
@assert all(in(keys(RATIO_OF_KM)), requested_km) "resolutions must be from $(sort(collect(keys(RATIO_OF_KM))))"

default_run && run_resolution(4//1; save_figure = false)   # discarded warm-up pass
results = [run_resolution(RATIO_OF_KM[km]) for km in requested_km]

#=
## Write `perf-scaling.csv` for `perf-scaling-lines.jl`
=#
const CSV_COLUMNS = (:resolution_km, :nx, :ny, :n_dof, :n_detached,
                     :elapsed_lin_cpu, :elapsed_lin_gpu, :elapsed_pt_cpu, :elapsed_pt_gpu,
                     :pt_cpu_iterations, :pt_gpu_iterations,
                     :pt_cpu_converged, :pt_gpu_converged)

#=
The default invocation truncates and writes the header; a resolution given on the command
line appends to whatever is already there, so the separate 4 km process lands in the same
file. Rows therefore arrive in invocation order — `perf-scaling-lines.jl` sorts by
`resolution_km` rather than assuming it.
=#
csv_path = joinpath(@__DIR__, "perf-scaling.csv")

#=
Merge rather than append: a resolution split across processes (4 km) contributes nothing
useful for the legs it did not run, and those must not overwrite a sibling process's real
numbers.

The merge is **per leg, not per value**. A skipped leg records `NaN` for its elapsed time but
`0` iterations and `false` converged — neither of which is `NaN`, so a naive value-wise merge
silently clobbers the other process's `480 / true` with `0 / false`. Each leg's elapsed time
therefore decides the fate of all of that leg's columns.
=#
const LEG_COLUMNS = (
    (:elapsed_pt_cpu,  (:elapsed_pt_cpu, :pt_cpu_iterations, :pt_cpu_converged)),
    (:elapsed_pt_gpu,  (:elapsed_pt_gpu, :pt_gpu_iterations, :pt_gpu_converged)),
    (:elapsed_lin_cpu, (:elapsed_lin_cpu,)),
    (:elapsed_lin_gpu, (:elapsed_lin_gpu,)),
)
colindex(c) = findfirst(==(c), CSV_COLUMNS)
didnotrun(v) = v isa AbstractFloat && isnan(v)

function merge_row(old, new)
    merged = copy(new)
    for (elapsed_col, cols) in LEG_COLUMNS
        if didnotrun(new[colindex(elapsed_col)])
            for c in cols
                merged[colindex(c)] = old[colindex(c)]
            end
        end
    end
    return merged
end

rows = Dict{Int,Vector{Any}}()
if !default_run && isfile(csv_path)
    prev = readdlm(csv_path, ',')
    for i in 2:size(prev, 1)
        rows[Int(prev[i, 1])] = collect(prev[i, :])
    end
end
for r in results
    km  = r.resolution_km
    new = Any[getproperty(r, c) for c in CSV_COLUMNS]
    rows[km] = haskey(rows, km) ? merge_row(rows[km], new) : new
end

open(csv_path, "w") do io
    writedlm(io, reshape(collect(string.(CSV_COLUMNS)), 1, :), ',')
    writedlm(io, permutedims(hcat((rows[km] for km in sort(collect(keys(rows)); rev = true))...)), ',')
end
println("\nwrote $csv_path")

@printf("\n%-6s %10s %10s %10s %10s %10s\n", "km", "DOFs", "lin CPU", "lin GPU", "PT CPU", "PT GPU")
for r in results
    @printf("%-6d %10d %9.2fs %9.2fs %9.2fs %9.2fs\n", r.resolution_km, r.n_dof,
           r.elapsed_lin_cpu, r.elapsed_lin_gpu, r.elapsed_pt_cpu, r.elapsed_pt_gpu)
end

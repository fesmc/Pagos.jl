#=

# Performance scaling: linear solve & Gershgorin-tuned PT vs. resolution

Wall-clock sweep of four DIVA momentum-solve configurations over 32/16/8 km:

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
through untouched. Coarsening is not a conservative remap: it turns `f_grnd` into a sub-grid
grounded *fraction* (harmless — `is_grounded` thresholds it at `> 0` anyway) and can dilute
thin or isolated ice below the `H_min = 0` `is_ice` threshold. At 32 km this already erases
all 48 detached iceberg cells the 8 km restart carries, so `n_detached` drops to 0 —
expected, not a bug.

`refine` (bilinear) is implemented and wired through `resample` for a future 4 km leg, but
**4 km is deliberately not in `RESOLUTIONS` yet**: it needs ~3.9 GiB of host *and* ~3.9 GiB
of device memory for its [`MechanicState`](@ref) alone, before the `lu` factors, which does
not fit comfortably on this machine. Note also that refining adds no information — a 4 km leg
would resolve nothing the 8 km field did not already contain, and bilinear interpolation
smears `H` across the ice margin — so it would be a size stress-test, not a 4 km simulation.

## Timing convention

`diva_update!`'s β_eff correction runs once per resolution *outside* every timed region,
matching `run_solve` (`helpers.jl`) and the linear-solve scripts. All four legs therefore time
only the core solve on an identical, already-DIVA-corrected problem: `A_active \\ b_active` /
CUDSS for the linear legs, `pseudo_transient!` for the PT legs. Each GPU leg runs a discarded
warm-up first (CUDSS symbolic analysis; PT kernel compilation past `SOLVER_KWARGS`' tuning
cadence), and the **whole 32 km resolution is run twice with the first pass discarded**, so
CPU-side JIT of the assembly loops and both solver paths is paid before anything is recorded.

!!! note "Memory"
    8 km is the heaviest leg currently run: ~1 GiB of host and ~1 GiB of device memory for its
    [`MechanicState`](@ref), plus the `lu` factors. All four resolutions run in one process, so
    each leg drops its state and forces a `GC.gc()` (and `CUDA.reclaim()`) before the next
    begins.
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

const SOLVER_KWARGS = (
    abstol = 1e-3, maxiter = 2000, ncheck = 20, printout_every = 50,
    pseudo_timestep = GershgorinPseudoTimeStep(cfl = 0.8),
    convergence = ScaledResidual(),
    friction_update = ActiveFrictionUpdate(),
    tuning = AutotunedDynamicRelaxation(cadence = 50),
)
const CST = Constants{Float32}()

depthavg_speed(mech) = Float32.(sqrt.(
    (@views (interior(mech.velocity.depthaverage_x)[1:(end - 1), :, 1] .+
             interior(mech.velocity.depthaverage_x)[2:end, :, 1]) ./ 2) .^ 2 .+
    (@views (interior(mech.velocity.depthaverage_y)[:, 1:(end - 1), 1] .+
             interior(mech.velocity.depthaverage_y)[:, 2:end, 1]) ./ 2) .^ 2))

#=
## One resolution: build the state, run all four solves, time each core solve
=#
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

    mech = MechanicState(grid)
    fill_from_grid!(mech.topography.thickness, H_ice)
    fill_from_grid!(mech.topography.surface, z_srf)
    fill_from_grid!(mech.material.viscosity_depthaveraged, visc_bar)
    fill_from_grid!(mech.friction.beta_eff, beta_eff)
    fill_from_grid!(mech.friction.beta, beta)
    fill_from_grid3d!(mech.material.viscosity, visc3d)
    setdata!(mech.velocity.depthaverage_x, zero(T))
    setdata!(mech.velocity.depthaverage_y, zero(T))

    solver_cpu = PseudoTransientSolver(grid; SOLVER_KWARGS...)
    diva_update!(mech, solver_cpu, rt, mask)   # β_eff correction — excluded from every timing

    #=
    ### Linear: coefficients (as in `linear-tuned-gpu.jl` / `res-sweep.jl`)
    =#
    Hf    = interior(mech.topography.thickness)[:, :, 1]
    sf    = interior(mech.topography.surface)[:, :, 1]
    viscf = interior(mech.material.viscosity_depthaveraged)[:, :, 1]
    befff = interior(mech.friction.beta_eff)[:, :, 1]

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

    function extract_speed(u_solved)
        lsd.u .= 0
        lsd.u[idx] .= T.(u_solved)
        ux = zeros(T, nx, ny); uy = zeros(T, nx, ny)
        velocity!(ux, uy, lsd)
        ux_c = similar(ux); uy_c = similar(uy)
        for i in 1:nx, j in 1:ny
            im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
            ux_c[i, j] = (ux[im1, j] + ux[i, j]) / 2
            uy_c[i, j] = (uy[i, jm1] + uy[i, j]) / 2
        end
        return T.(sqrt.(ux_c .^ 2 .+ uy_c .^ 2))
    end

    #=
    ### 1. Linear, CPU
    =#
    elapsed_lin_cpu = @elapsed u_cpu = A_active \ b_active
    speed_lin_cpu = on_ice(extract_speed(u_cpu))
    @printf("  linear CPU: %8.2f s  (%d active DOFs)\n", elapsed_lin_cpu, length(idx))

    #=
    ### 2. Linear, GPU (CUDSS)

    `lsd_gpu` is dropped and the device pool reclaimed after each call, so the CUDSS factors
    are gone before the PT leg's `mech_gpu` is allocated rather than sharing the card with it.
    =#
    run_gpu = CUDA.functional()
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
    ### 3. Gershgorin PT, CPU — reuses `mech` (β_eff corrected, velocity zeroed above)
    =#
    elapsed_pt_cpu = @elapsed res_pt_cpu =
        pseudo_transient!(mech, CST, solver_cpu, rt, DIVAMomentumBalance(), mask)
    speed_pt_cpu = on_ice(depthavg_speed(mech))
    @printf("  PT CPU:     %8.2f s  (%d it, converged = %s)\n",
           elapsed_pt_cpu, res_pt_cpu.iterations, res_pt_cpu.converged)

    #=
    ### 4. Gershgorin PT, GPU
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
        elapsed_pt_gpu = NaN
        res_pt_gpu = (; converged = false, iterations = 0)
        speed_pt_gpu = fill(NaN32, nx, ny)
    end

    save_figure && save_resolution_figure(; resolution_km, ratio, nx, ny, z_bed,
        elapsed_lin_cpu, elapsed_lin_gpu, elapsed_pt_cpu, elapsed_pt_gpu,
        res_pt_cpu, res_pt_gpu, run_gpu,
        speed_lin_cpu, speed_lin_gpu, speed_pt_cpu, speed_pt_gpu)

    out = (; resolution_km, nx, ny, n_dof = length(idx), n_detached,
          elapsed_lin_cpu, elapsed_lin_gpu, elapsed_pt_cpu, elapsed_pt_gpu,
          pt_cpu_iterations = res_pt_cpu.iterations, pt_gpu_iterations = res_pt_gpu.iterations,
          pt_cpu_converged = res_pt_cpu.converged, pt_gpu_converged = res_pt_gpu.converged)

    mech = nothing; lsd = nothing; masks = nothing
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

Add `1//2` to extend the sweep to 4 km — see the memory note in the module docstring first.
=#
const RESOLUTIONS = [4//1, 2//1, 1//1]   # → 32, 16, 8 km

run_resolution(4//1; save_figure = false)      # discarded warm-up pass
results = [run_resolution(r) for r in RESOLUTIONS]

#=
## Write `perf-scaling.csv` for `perf-scaling-lines.jl`
=#
const CSV_COLUMNS = (:resolution_km, :nx, :ny, :n_dof, :n_detached,
                     :elapsed_lin_cpu, :elapsed_lin_gpu, :elapsed_pt_cpu, :elapsed_pt_gpu,
                     :pt_cpu_iterations, :pt_gpu_iterations,
                     :pt_cpu_converged, :pt_gpu_converged)

csv_path = joinpath(@__DIR__, "perf-scaling.csv")
open(csv_path, "w") do io
    writedlm(io, reshape(collect(string.(CSV_COLUMNS)), 1, :), ',')
    writedlm(io, [getproperty(r, c) for r in results, c in CSV_COLUMNS], ',')
end
println("\nwrote $csv_path")

@printf("\n%-6s %10s %10s %10s %10s %10s\n", "km", "DOFs", "lin CPU", "lin GPU", "PT CPU", "PT GPU")
for r in results
    @printf("%-6d %10d %9.2fs %9.2fs %9.2fs %9.2fs\n", r.resolution_km, r.n_dof,
           r.elapsed_lin_cpu, r.elapsed_lin_gpu, r.elapsed_pt_cpu, r.elapsed_pt_gpu)
end

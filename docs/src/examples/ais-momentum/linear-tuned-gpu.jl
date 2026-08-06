#=

# DIVA: direct linear solve vs. fixed-tune PT vs. Gershgorin-tuned PT vs. GPU

Four solves of the same DIVA problem, all on the 8 km restart:

 1. **Linear** — a direct sparse solve ([`LinearMomentumSolver2D`](@ref)), CPU, Float32.
 2. **Fixed-tune PT** — [`PseudoTransientSolver`](@ref) with [`ViscosityPseudoTimeStep`](@ref)
    (Sandip et al. 2024 Eq. 7, `Δτ ∝ 1/η`), CPU, Float32.
 3. **Gershgorin-tuned PT** — the library default, [`GershgorinPseudoTimeStep`](@ref), CPU,
    Float32.
 4. **Gershgorin-tuned PT** — same as 3, on GPU.

3 vs. 4 is the same comparison `cpu-gpu.jl` makes (there at Float64); this script repeats it
at Float32 as one leg of a wider comparison rather than as its own subject. 2 vs. 3 compares
`pseudo_timestep` (`Δτ` rule) *and* `tuning` (relaxation strategy) together, since
`AutotunedDynamicRelaxation` requires [`GershgorinPseudoTimeStep`](@ref) and cannot be held
fixed across both legs — so this ends up close to `hand-auto-tuned.jl`'s `FixedTuning()` vs.
`AutotunedDynamicRelaxation()` comparison, with `pseudo_timestep` riding along as a second,
inseparable axis of difference rather than an isolated one.

## The direct linear solve: what it means here, and why the mask matters

[`LinearMomentumSolver2D`](@ref) assembles the same DIVA operator PT iterates toward and
factorizes it directly (`lu`), rather than relaxing to it — see
`benchmark/basics/pseudotransient/slab-pt-vs-linear.jl` for the same comparison on a uniform
slab. It only ships a `RegularGrid`/periodic-boundary constructor and knows nothing about
[`momentum_mask!`](@ref)'s `is_momentum_solved`, so three things below are new code that the
constructor call itself does not give for free:

 1. **Coefficients.** `N` (`= 2 H μ̄`, at cell centers), `N_ab` (the same quantity harmonically
    interpolated in viscosity / arithmetically in `H` onto the corners), `β_acx`/`β_acy`
    (arithmetic face-average of the DIVA-corrected `β_eff` from [`diva_update!`](@ref)) and
    `τ_x`/`τ_y` (`ρgH_face·∂s/∂x`, matching [`drivingstress!`](@ref)'s sign convention) are
    built by hand from the same restart fields `helpers.jl` already loaded, on a Float32
    scratch [`MechanicState`](@ref) whose only job is to run [`diva_update!`](@ref) once.
    Cross-checked two independent ways — this arithmetic, and calling Pagos' own `lerp`/`hlerp`/
    [`drivingstress!`](@ref) on the Chmy-native fields with the one-cell index shift their
    "vertex `i` sits between centres `i-1` and `i`" convention needs to line up with this
    solver's "face `i` sits between centres `i` and `i+1`" one — and against
    [`GershgorinPseudoTimeStep`](@ref)'s own row-sum formula, which re-derives the same `N`/`N_ab`
    coefficients independently. All three agree.
 2. **The mask.** The full periodic system is singular over most of the domain: open ocean has
    no friction and no driving stress (`β = 0`, `τ_d = 0`, pure diffusion with no source), and
    detached icebergs have no basal drag and no membrane connection to anything — see
    [`momentum_mask!`](@ref)'s docstring for why that is unsolvable, not merely uninteresting.
    The DOFs are restricted to `is_momentum_solved` before factorizing, requiring **both**
    flanking cells of a face to be solved (the strict, `node_fully_active`-style rule — the
    permissive "either" rule leaves faces whose row is built almost entirely from an excluded
    neighbour, which are close enough to singular that `lu` still returns a "solution" dominated
    by round-off). Excluded DOFs stay at their initial `0`, the same fallback PT uses.
 3. **Extraction.** `LinearMomentumSolver2D`'s own `velocity!(ux, uy, lsd)` unpacks the solved
    vector; a small local average converts the `(nx,ny)` face-valued `ux`/`uy` (face `i` between
    centres `i`, `i+1`) to a cell-centred speed, matching the on-ice masking every other script
    here applies before plotting.

!!! warning "The direct solve does not reproduce the PT solution on this geometry"
    Validated on a small synthetic DIVA problem with smoothly-varying `H`, viscosity *and*
    friction (correlation 0.9997, mean ratio 0.998 against a Chmy-native PT solve of the same
    fields) — the coefficient construction and mask restriction are correct in that setting.
    On the real 8 km restart, restricted to the identical `is_momentum_solved` DOFs, it diverges
    substantially from the converged Gershgorin-tuned PT solution: not a handful of outlier
    cells (an ice-thickness floor on top of the mask, up to 200 m, moves the discrepancy
    negligibly) but a pervasive one, ~195% median relative difference even after trimming the
    worst 5% by absolute error. The sparse factorization's own residual is essentially zero
    throughout — which only means self-consistent, not correct: a tiny `lu` residual is normal
    for a badly-conditioned system and is no evidence of an accurate solution. Two explanations
    remain open and are **not** distinguished here: a genuine conditioning problem from directly
    inverting DIVA across real grounding lines (`β` can jump by orders of magnitude between
    adjacent cells, unlike the smooth synthetic check) and complex domain connectivity, or a
    residual bug in `LinearMomentumSolver2D`'s legacy `loop1_coeffs`/`loop2_coeffs` assembly that
    the uniform-slab test (`test/mechanics/slab.jl`) cannot expose, since most of its terms are
    insensitive to spatially-varying coefficients on a uniform field. Read the "Linear" panel
    below as a demonstration of *that gap*, not as a working direct-solve baseline — it is why
    [`GershgorinPseudoTimeStep`](@ref) is the library default for real geometry rather than a
    convenience.
=#
resolution_km = 8

using Pkg
Pkg.activate(joinpath(@__DIR__, "../../.."))
using CUDA
if CUDA.functional()
    CUDA.zeros(1)
    CUDA.synchronize()
end

include(joinpath(@__DIR__, "helpers.jl"))
using SparseArrays, LinearAlgebra

#=
## 1. Linear: masked direct solve, CPU, Float32

`grid32`/`rt32` are reused below for the two CPU PT legs too, since `diva_update!`'s
`β_eff` is the only thing that needs a live [`MechanicState`](@ref) here; `run_solve` builds
and tears down its own.
=#
grid32 = StaggeredGrid(Float32, lx, ly, dx, dy, layering)
rt32   = Runtime(grid32)

mech32 = MechanicState(grid32)
fill_from_grid!(mech32.topography.thickness, Float32.(H_ice))
fill_from_grid!(mech32.topography.surface, Float32.(z_srf))
fill_from_grid!(mech32.material.viscosity_depthaveraged, Float32.(visc_bar))
fill_from_grid!(mech32.friction.beta_eff, Float32.(beta_eff))
fill_from_grid!(mech32.friction.beta, Float32.(beta))
fill_from_grid3d!(mech32.material.viscosity, Float32.(visc3d))
setdata!(mech32.velocity.depthaverage_x, 0.0f0)
setdata!(mech32.velocity.depthaverage_y, 0.0f0)

solver32 = PseudoTransientSolver(grid32; SOLVER_KWARGS...)
diva_update!(mech32, solver32, rt32, mask)   # writes the DIVA-corrected β_eff

H32    = interior(mech32.topography.thickness)[:, :, 1]
s32    = interior(mech32.topography.surface)[:, :, 1]
visc32 = interior(mech32.material.viscosity_depthaveraged)[:, :, 1]
beff32 = interior(mech32.friction.beta_eff)[:, :, 1]
mech32 = nothing
GC.gc()

i_idx = PeriodicIndexing(1, nx)
j_idx = PeriodicIndexing(1, ny)
rho_ice = Float32(cst.density_ice)
grav    = Float32(cst.gravity)

N     = similar(H32); N_ab = similar(H32)
b_acx = similar(H32); b_acy = similar(H32)
tx    = similar(H32); ty   = similar(H32)
for i in 1:nx, j in 1:ny
    im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
    N[i, j] = 2 * H32[i, j] * visc32[i, j]
    eta_ab = 4 / (inv(visc32[i, j]) + inv(visc32[ip1, j]) + inv(visc32[i, jp1]) + inv(visc32[ip1, jp1]))
    H_ab   = (H32[i, j] + H32[ip1, j] + H32[i, jp1] + H32[ip1, jp1]) / 4
    N_ab[i, j] = 2 * eta_ab * H_ab
    b_acx[i, j] = (beff32[i, j] + beff32[ip1, j]) / 2
    b_acy[i, j] = (beff32[i, j] + beff32[i, jp1]) / 2
    Hx = (H32[i, j] + H32[ip1, j]) / 2
    Hy = (H32[i, j] + H32[i, jp1]) / 2
    tx[i, j] = rho_ice * grav * Hx * (s32[ip1, j] - s32[i, j]) / dx
    ty[i, j] = rho_ice * grav * Hy * (s32[i, jp1] - s32[i, j]) / dy
end

# `(nx-1)*dx`, not `nx*dx`: `RegularGrid` builds its axis as `0:dx:lx`, so this is what
# gives exactly `nx` grid points (an off-by-one otherwise) -- see `test/mechanics/slab.jl`.
regular_grid = RegularGrid(Float32, (nx - 1) * dx, (ny - 1) * dy, dx, dy)
lsd = LinearMomentumSolver2D(regular_grid, DIVAMomentumBalance())
ux0 = zeros(Float32, nx, ny); uy0 = zeros(Float32, nx, ny)
populate_vectors!(lsd, N, N_ab, ux0, uy0, tx, ty, b_acx, b_acy)

is_solved = interior(topo.mask.is_momentum_solved)[:, :, 1]
active_dof = falses(2 * nx * ny)
for i in 1:nx, j in 1:ny
    im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
    active_dof[Pagos._ij2n_ux(i, j, nx, ny)] = is_solved[i, j] & is_solved[ip1, j]
    active_dof[Pagos._ij2n_uy(i, j, nx, ny)] = is_solved[i, j] & is_solved[i, jp1]
end
idx = findall(active_dof)
@printf("Linear: %d / %d DOFs active (is_momentum_solved, both flanking cells)\n",
       length(idx), 2 * nx * ny)

A_active = lsd.A[idx, idx]
b_active = lsd.b[idx]
elapsed_lin = @elapsed u_active = A_active \ b_active
lsd.u .= 0
lsd.u[idx] .= Float32.(u_active)
ux_lin = zeros(Float32, nx, ny); uy_lin = zeros(Float32, nx, ny)
velocity!(ux_lin, uy_lin, lsd)

ux_c = similar(ux_lin); uy_c = similar(uy_lin)
for i in 1:nx, j in 1:ny
    im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
    ux_c[i, j] = (ux_lin[im1, j] + ux_lin[i, j]) / 2
    uy_c[i, j] = (uy_lin[i, jm1] + uy_lin[i, j]) / 2
end
speed_lin = Float32.(sqrt.(ux_c .^ 2 .+ uy_c .^ 2))
println("Linear, CPU Float32: elapsed = ", round(elapsed_lin, digits = 2), "s ",
       "(note: SparseArrays' `lu` promotes Float32 sparse matrices to Float64 internally --",
       " a SuiteSparse/UMFPACK limitation, not a Pagos one)")

#=
## 2. Fixed-tune PT (`ViscosityPseudoTimeStep`) vs. 3. Gershgorin-tuned PT — CPU, Float32

`AutotunedDynamicRelaxation` (`SOLVER_KWARGS`'s own default `tuning`) requires
[`GershgorinPseudoTimeStep`](@ref) (`_check_tuning` throws otherwise), so the fixed-tune leg
below also swaps `tuning` to `FixedTuning()` — see the module docstring for why that makes
this a two-preset comparison rather than a single isolated variable.

`ViscosityPseudoTimeStep` ignores basal drag when choosing `Δτ` (its own docstring's warning),
which under real, friction-dominated grounded ice stalls convergence rather than merely
slowing it: on this restart the residual plateaus around `0.93` (abstol is `1e-3`) within the
first 20 iterations and stays there, so `maxiter` is capped at a few hundred rather than chased
up toward convergence that inspection shows will not arrive. The Gershgorin-tuned legs below
use `SOLVER_KWARGS` unchanged, i.e. the library default (`AutotunedDynamicRelaxation(cadence =
50)`), same as every other script in this folder.
=#
diva_fixed = run_solve(DIVAMomentumBalance(), grid32, rt32, mask;
    (; SOLVER_KWARGS..., pseudo_timestep = ViscosityPseudoTimeStep(Float32),
      tuning = FixedTuning(), maxiter = 500)...)
println("Fixed-tune PT, CPU Float32: ",
       (; diva_fixed.converged, diva_fixed.iterations, diva_fixed.elapsed, diva_fixed.residual))

diva_gersh = run_solve(DIVAMomentumBalance(), grid32, rt32, mask; SOLVER_KWARGS...)
println("Gershgorin-tuned PT, CPU Float32: ",
       (; diva_gersh.converged, diva_gersh.iterations, diva_gersh.elapsed, diva_gersh.residual))

#=
## 4. Gershgorin-tuned PT — GPU, Float32

Same structure as `cpu-gpu.jl` (CUDA touched before `helpers.jl`'s `set_theme!` call, a
discarded warm-up solve run in place on `mech_gpu` past `SOLVER_KWARGS`'s tuning cadence) —
see that script for the full rationale. `SOLVER_KWARGS` unchanged here too, matching the CPU
leg above.
=#
if CUDA.functional()
    arch_gpu = Arch(CUDABackend())
    grid_gpu = StaggeredGrid(arch_gpu, Float32, lx, ly, dx, dy, layering)
    rt_gpu   = Runtime(grid_gpu)

    mech_cpu_for_gpu = MechanicState(grid32)
    fill_from_grid!(mech_cpu_for_gpu.topography.thickness, Float32.(H_ice))
    fill_from_grid!(mech_cpu_for_gpu.topography.surface, Float32.(z_srf))
    fill_from_grid!(mech_cpu_for_gpu.material.viscosity_depthaveraged, Float32.(visc_bar))
    fill_from_grid!(mech_cpu_for_gpu.friction.beta_eff, Float32.(beta_eff))
    fill_from_grid!(mech_cpu_for_gpu.friction.beta, Float32.(beta))
    fill_from_grid3d!(mech_cpu_for_gpu.material.viscosity, Float32.(visc3d))
    setdata!(mech_cpu_for_gpu.velocity.depthaverage_x, 0.0f0)
    setdata!(mech_cpu_for_gpu.velocity.depthaverage_y, 0.0f0)

    topo_gpu = Pagos.Adapt.adapt(CuArray, topo)
    mask_gpu = IceMask(topo_gpu.mask.is_momentum_solved)

    function solve_on!(mech, grid, rt, mask; solver_kwargs...)
        solver = PseudoTransientSolver(grid; solver_kwargs...)
        diva_update!(mech, solver, rt, mask)
        elapsed = @elapsed result = pseudo_transient!(mech, Constants{Float32}(), solver, rt,
                                                       DIVAMomentumBalance(), mask)
        speed = Array(Float32.(sqrt.(
            (@views (interior(mech.velocity.depthaverage_x)[1:(end - 1), :, 1] .+
                     interior(mech.velocity.depthaverage_x)[2:end, :, 1]) ./ 2) .^ 2 .+
            (@views (interior(mech.velocity.depthaverage_y)[:, 1:(end - 1), 1] .+
                     interior(mech.velocity.depthaverage_y)[:, 2:end, 1]) ./ 2) .^ 2)))
        return (; speed, elapsed, iterations = result.iterations, converged = result.converged,
               residual = result.residual)
    end

    mech_gpu = Pagos.Adapt.adapt(CuArray, mech_cpu_for_gpu)
    mech_cpu_for_gpu = nothing
    GC.gc()
    CUDA.reclaim()

    solve_on!(mech_gpu, grid_gpu, rt_gpu, mask_gpu;
        (; SOLVER_KWARGS..., maxiter = 55, abstol = 0.0)...)
    setdata!(mech_gpu.velocity.depthaverage_x, 0.0f0)
    setdata!(mech_gpu.velocity.depthaverage_y, 0.0f0)
    println("GPU warm-up done (kernels compiled, not timed).")

    diva_gpu = solve_on!(mech_gpu, grid_gpu, rt_gpu, mask_gpu; SOLVER_KWARGS...)
    mech_gpu = nothing
    GC.gc()
    println("Gershgorin-tuned PT, GPU Float32: ",
           (; diva_gpu.converged, diva_gpu.iterations, diva_gpu.elapsed, diva_gpu.residual))
else
    println("CUDA.functional() == false -- GPU panel will be blank.")
    diva_gpu = (; speed = fill(NaN32, nx, ny), elapsed = NaN, iterations = 0,
               converged = false, residual = NaN)
end

#=
## Diagnostics

`speed_lin` vs. `diva_gersh.speed` is the divergence flagged in the warning above, not a
useful accuracy statement about either solve; `diva_gersh.speed` vs. `diva_gpu.speed` is the
CPU/GPU consistency check `cpu-gpu.jl` already makes, repeated here at Float32.
=#
function pairwise_stats(name, a, b)
    d = on_ice(a .- b)
    vals = filter(!isnan, d)
    @printf("%s (on-ice, %d cells): RMSE = %.4g m/yr, median|Δ| = %.4g m/yr, max|Δ| = %.4g m/yr\n",
           name, length(vals), sqrt(mean(abs2, vals)), median(abs.(vals)), maximum(abs, vals))
    return d
end

diff_lin_gersh = pairwise_stats("Linear vs. Gershgorin PT (CPU)", speed_lin, diva_gersh.speed)
diff_fixed_gersh = pairwise_stats("Fixed-tune vs. Gershgorin PT (CPU)", diva_fixed.speed, diva_gersh.speed)
diff_cpu_gpu = CUDA.functional() ?
    pairwise_stats("Gershgorin PT, CPU vs. GPU", diva_gersh.speed, diva_gpu.speed) :
    fill(NaN32, nx, ny)

#=
## Figure

MEaSUREs-style velocity colouring (white → dodgerblue → yellow → red → darkred, log-scaled,
each colour landing exactly on 0/100/400/700/1000 m/yr) over a `:bukavu` bathymetry backdrop
drawn first in each panel, so the ocean/off-ice area (transparent `NaN` in the speed layer)
reads as seafloor rather than blank. One row, four columns: Linear, then the three PT
variants. Cropped by 30 cells in `x` and 80 in `y` (each side) to trim the mostly-empty
domain margin around the ice sheet.
=#
crop_x, crop_y = 30, 80
ix = (crop_x + 1):(nx - crop_x)
iy = (crop_y + 1):(ny - crop_y)
xc_c, yc_c = xc[ix], yc[iy]
crop(data) = data[ix, iy]

stops = [0, 20, 100, 400, 700, 1000]
speed_cmap = cgrad(
    [:white, :white, :dodgerblue4, :lightgoldenrod1, :orangered, :darkred],
    stops ./ 1000,
)
z_bed_c = crop(Float32.(z_bed))
crange = extrema(stops)

bmap = cgrad([ :white])
perf_label(r) = @sprintf("%d it · %.1f s%s", r.iterations, r.elapsed,
                         r.converged ? "" : " (not converged)")

function panel!(fig_pos, title, subtitle, data)
    ax = Axis(fig_pos, title = title, subtitle = subtitle, aspect = DataAspect())
    heatmap!(ax, xc_c, yc_c, z_bed_c; colormap = :oleron, colorrange = (-6000, 6000))
    hm = heatmap!(ax, xc_c, yc_c, max.(crop(data), 1); colorrange = crange,
        colormap = speed_cmap,
        lowclip = speed_cmap[1], highclip = speed_cmap[end])
    hidedecorations!(ax)
    return hm
end

fig = Figure(size = (1650, 500), fontsize = 18)

panel!(fig[1, 1], "Linear (masked direct solve)",
    @sprintf("%.1f s (diverges from PT, see warning)", elapsed_lin), on_ice(speed_lin))
panel!(fig[1, 2], "Fixed-tune PT", perf_label(diva_fixed), on_ice(diva_fixed.speed))
panel!(fig[1, 3], "Gershgorin-tuned PT, CPU", perf_label(diva_gersh), on_ice(diva_gersh.speed))
hm_gpu = panel!(fig[1, 4], "Gershgorin-tuned PT, GPU", perf_label(diva_gpu), on_ice(diva_gpu.speed))

Colorbar(fig[2, 2:3], hm_gpu, vertical = false, width = Relative(0.5), flipaxis = false, height = Relative(1), label = "speed (m/yr)")
rowsize!(fig.layout, 2, 20)
colgap!(fig.layout, 5)
save("$figdir/linear-tuned-gpu-$(resolution_km)km.png", fig)
fig

#=

# DIVA: masked direct linear solve, CPU vs. GPU (CUDSS) — first try, 8 km

A first try of the *externalised* GPU-compatible linear solve on real AIS geometry:
[`LinearMomentumSolver2D`](@ref)'s `velocity!` dispatches to a CUDSS-backed direct solve
(`ext/PagosCUDSSExt.jl`) whenever its buffers are `CuVector`/`CuSparseMatrixCSR` — the same
struct, kernels ([`populate_vectors!`](@ref)) and public API as the CPU path, no separate
GPU code to maintain. `benchmark/basics/pseudotransient/slab-pt-vs-linear.jl` exercises this
on a synthetic uniform slab; this script is the first attempt at the real, masked AIS
problem, at the resolution `linear-tuned-gpu.jl` already validated on CPU (8 km). The name
`linsolve-gpu.jl` is aspirational: once this leg is trusted, the plan is to loop it (and the PT
legs) over the restarts at other resolutions — not attempted yet.

## What's reused vs. new here

Coefficient assembly (`N`, `N_ab`, `β_acx`/`β_acy`, `τ_x`/`τ_y`) and the DOF mask
(`is_momentum_solved`, permissive either-flanking-cell rule) are copied verbatim from
`linear-tuned-gpu.jl` — see that script's module docstring for the full derivation, the
validation against Pagos' own residual kernels, and the two-bug history (`loop1_coeffs` sign,
`N` factor-of-2) behind why the numbers below are trustworthy. The only new part is the
solve itself: instead of `A_active \\ b_active` (CPU `SparseArrays.lu`), the reduced system
moves to the GPU as `CuSparseMatrixCSR`/`CuVector` and solves through
[`LinearMomentumSolver2D`](@ref)'s own `velocity!` — exercising the *same* code path
`slab-pt-vs-linear.jl` benchmarks, just with `A`/`u`/`b` built by hand rather than by
`populate_vectors!`.

!!! note "Environment: `CUDSS` compat had to widen"
    `CUDSS` was a weak dependency already (`ext/PagosCUDSSExt.jl`), but unreachable in
    practice: `Project.toml` pinned `CUDSS = "0.7"`, and CUDSS 0.7 depends on `CUDACore`,
    which only pairs with CUDA.jl ≥ 6 — one major version past the `CUDA = "5"` this package
    (and `Chmy`) is pinned to. CUDSS 0.6.5 depends on `CUDA = "5.4.0"` directly (no
    `CUDACore` split yet) and exposes the identical `CudssSolver`/`cudss` API
    `PagosCUDSSExt.jl` calls, so the fix was widening the bound to `CUDSS = "0.6 - 0.7"`
    rather than touching `CUDA`'s own pin — no code changes needed, confirmed by reproducing
    a tiny hand-built `CuSparseMatrixCSR` solve against the CPU answer to Float64 round-off
    before running this script for real.

!!! warning "The DOF reduction itself is still CPU-only"
    `A_active = lsd.A[idx, idx]` is CPU sparse fancy-indexing — cheap here (8 km, ~1.16M DOFs
    before reduction), but it means this script moves an *already-reduced* system to the GPU
    rather than assembling the reduced system there directly. A real GPU-resident pipeline
    would need either a GPU-native reduction or a regularization scheme that keeps the full
    periodic system's shape (e.g. an identity row/RHS-zero on excluded DOFs) so
    `populate_vectors!`'s existing kernels could write the whole system without a CPU round
    trip. Left open, same as the `N_ab`-masking gap `linear-tuned-gpu.jl` flags.
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
using SparseArrays, LinearAlgebra, CUDA.CUSPARSE, CUDSS

#=
## Coefficients and the masked DOF system (copied from `linear-tuned-gpu.jl`)
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
    N[i, j] = H32[i, j] * visc32[i, j]
    eta_ab = 4 / (inv(visc32[i, j]) + inv(visc32[ip1, j]) + inv(visc32[i, jp1]) + inv(visc32[ip1, jp1]))
    H_ab   = (H32[i, j] + H32[ip1, j] + H32[i, jp1] + H32[ip1, jp1]) / 4
    N_ab[i, j] = eta_ab * H_ab
    b_acx[i, j] = (beff32[i, j] + beff32[ip1, j]) / 2
    b_acy[i, j] = (beff32[i, j] + beff32[i, jp1]) / 2
    Hx = (H32[i, j] + H32[ip1, j]) / 2
    Hy = (H32[i, j] + H32[i, jp1]) / 2
    tx[i, j] = rho_ice * grav * Hx * (s32[ip1, j] - s32[i, j]) / dx
    ty[i, j] = rho_ice * grav * Hy * (s32[i, jp1] - s32[i, j]) / dy
end

regular_grid = RegularGrid(Float32, (nx - 1) * dx, (ny - 1) * dy, dx, dy)
lsd = LinearMomentumSolver2D(regular_grid, DIVAMomentumBalance())
ux0 = zeros(Float32, nx, ny); uy0 = zeros(Float32, nx, ny)
populate_vectors!(lsd, N, N_ab, ux0, uy0, tx, ty, b_acx, b_acy)

is_solved = interior(masks.is_momentum_solved)[:, :, 1]
active_dof = falses(2 * nx * ny)
for i in 1:nx, j in 1:ny
    im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
    active_dof[Pagos._ij2n_ux(i, j, nx, ny)] = is_solved[i, j] | is_solved[ip1, j]
    active_dof[Pagos._ij2n_uy(i, j, nx, ny)] = is_solved[i, j] | is_solved[i, jp1]
end
idx = findall(active_dof)
@printf("%d / %d DOFs active (is_momentum_solved, either flanking cell)\n",
       length(idx), 2 * nx * ny)

A_active = lsd.A[idx, idx]
b_active = lsd.b[idx]

function extract_speed(u_solved)
    lsd.u .= 0
    lsd.u[idx] .= Float32.(u_solved)
    ux = zeros(Float32, nx, ny); uy = zeros(Float32, nx, ny)
    velocity!(ux, uy, lsd)
    ux_c = similar(ux); uy_c = similar(uy)
    for i in 1:nx, j in 1:ny
        im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
        ux_c[i, j] = (ux[im1, j] + ux[i, j]) / 2
        uy_c[i, j] = (uy[i, jm1] + uy[i, j]) / 2
    end
    return Float32.(sqrt.(ux_c .^ 2 .+ uy_c .^ 2))
end

#=
## CPU solve (baseline, matches `linear-tuned-gpu.jl`)
=#
elapsed_cpu = @elapsed u_cpu = A_active \ b_active
println("Linear, CPU Float32 (SparseArrays lu): elapsed = ", round(elapsed_cpu, digits = 2), "s")

#=
## GPU solve, via the externalised CUDSS-backed `LinearMomentumSolver2D`

`populate_vectors!` is not called here — `A_active`/`b_active` are already fully assembled
(reduced, CPU-side, see the warning above) — so this only needs the struct fields
`velocity!`'s CUDSS dispatch actually reads: `.A`, `.u`, `.b`, `.solver_cache`. The rest
(`dynamics`, `nx`, `ny`, the `dxdx_`/`dydy_`/`dxdy_` stencil scalars, `perm`, `i_idx`,
`j_idx`) exist for [`populate_vectors!`](@ref)'s full-grid assembly, which this call never
reaches, so they're placeholders.
=#
function gpu_cudss_solve(A_active, b_active)
    T = eltype(b_active)
    A_gpu = CuSparseMatrixCSR(A_active)
    b_gpu = CuArray(b_active)
    u_gpu = CUDA.zeros(T, length(b_active))
    lsd_gpu = LinearMomentumSolver2D(
        DIVAMomentumBalance(), 0, 0, zero(T), zero(T), zero(T),
        u_gpu, similar(u_gpu), b_gpu, A_gpu, CUDA.zeros(Int, 1),
        PeriodicIndexing(1, 1), PeriodicIndexing(1, 1), Ref{Any}(nothing),
    )
    velocity!(lsd_gpu)   # PagosCUDSSExt: analysis (cached) + factorization + solve
    CUDA.synchronize()
    return Array(lsd_gpu.u)
end

if CUDA.functional()
    gpu_cudss_solve(A_active, b_active)   # warm-up: CUDSS analysis + kernel compile
    elapsed_gpu = @elapsed u_gpu = gpu_cudss_solve(A_active, b_active)
    println("Linear, GPU Float32 (CUDSS lu):        elapsed = ", round(elapsed_gpu, digits = 2), "s")

    diff = abs.(Float64.(u_cpu) .- Float64.(u_gpu))
    @printf("CPU vs. GPU linear solve (%d active DOFs): max|Δ| = %.4g, median|Δ| = %.4g\n",
           length(idx), maximum(diff), median(diff))
else
    println("CUDA.functional() == false -- GPU leg skipped.")
    u_gpu = copy(u_cpu)
    elapsed_gpu = NaN
end

#=
## Figure

Same panel/colorbar layout as `linear-tuned-gpu.jl`: one row, two columns, colorbar centred
below both.
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

function panel!(fig_pos, title, subtitle, data)
    ax = Axis(fig_pos, title = title, subtitle = subtitle, aspect = DataAspect())
    heatmap!(ax, xc_c, yc_c, z_bed_c; colormap = :oleron, colorrange = (-6000, 6000))
    hm = heatmap!(ax, xc_c, yc_c, max.(crop(data), 1); colorrange = crange,
        colormap = speed_cmap,
        lowclip = speed_cmap[1], highclip = speed_cmap[end])
    hidedecorations!(ax)
    return hm
end

fig = Figure(size = (900, 500), fontsize = 18)

panel!(fig[1, 1], "Linear (CPU, SparseArrays lu)",
    @sprintf("%.1f s", elapsed_cpu), on_ice(extract_speed(u_cpu)))
hm_gpu = panel!(fig[1, 2], "Linear (GPU, CUDSS lu)",
    CUDA.functional() ? @sprintf("%.2f s", elapsed_gpu) : "CUDA unavailable",
    on_ice(extract_speed(u_gpu)))

Colorbar(fig[2, 1:2], hm_gpu, vertical = false, width = Relative(0.5), flipaxis = false, height = Relative(1), label = "speed (m/yr)")
rowsize!(fig.layout, 2, 20)
colgap!(fig.layout, 5)
save("$figdir/linsolve-gpu-$(resolution_km)km.png", fig)
fig

#=

# DIVA: Pagos vs. Yelmo — surface speed, RMSE only

!!! note "Vertical layering: index-matched to Yelmo, not value-matched"
    Yelmo's `zeta`/`zeta_ac` do not correspond to [`QuadraticSigmaTransform`](@ref)'s own
    derivation (checked numerically — different values at every interior level), and
    reproducing them exactly would mean bypassing the parametric transform for a
    `CorrectedVerticalLayering` built straight from Yelmo's arrays. Not worth it here: this
    script uses a plain `QuadraticSigmaTransform` grid with Yelmo's own layer *count*
    (`nz = 11`) and copies Yelmo's per-layer viscosity in by index, not by matching sigma
    value. That is why this comparison is surface-only (a quantity insensitive to the exact
    interior discretization) rather than a profile comparison.

`ux`/`uy` are read here, not in `helpers.jl`'s loading section, as a **2D surface slice
straight from the file** (`ds["ux"][:, :, end, 1]`) — `zeta[end] == 1.0` is the surface (the
`zeta` printed while writing this script runs bed-to-surface,
`[0.0, 0.01, ..., 1.0]`). The other 10 layers of the full 3D `ux`/`uy` are never read into
memory at all, which is the whole point: this comparison, per its own scoping decision, is
surface-only, so there is nothing to gain from the other layers and no reason to pay for
them. Float16 is safe here — Yelmo's surface speeds top out under 4300 m/yr in this file,
comfortably inside Float16's ~65504 magnitude cap (checked, like every other narrow-type
load in `helpers.jl`, against the file's own values rather than assumed).
=#
include(joinpath(@__DIR__, "helpers.jl"))

ux_srf, uy_srf = NCDataset(restart_file) do ds
    (Float16.(ds["ux"][:, :, end, 1]), Float16.(ds["uy"][:, :, end, 1]))
end
speed_yelmo_srf = on_ice(sqrt.(Float32.(ux_srf) .^ 2 .+ Float32.(uy_srf) .^ 2))

# The depth-averaged speed [`run_solve`](@ref) returns is not the surface one — comparing it
# against Yelmo's surface speed would be comparing two different physical quantities.
# `velocities3D!`'s Eq. 17 surface reconstruction is what belongs on this side of the
# comparison, so this section solves DIVA directly (rather than through `run_solve`) and
# reconstructs the surface velocity from the converged state before comparing.
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
speed_pagos_srf = on_ice(Float32.(sqrt.(
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
save("$figdir/yelmo-pagos.png", fig5)
fig5

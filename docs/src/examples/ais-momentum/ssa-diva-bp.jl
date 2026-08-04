#=

# SSA vs. DIVA vs. Blatter-Pattyn

!!! warning "Not quite apples-to-apples on friction"
    The SSA run uses Yelmo's own `beta_eff` directly but that field
    is *already* Yelmo's own DIVA vertical-shear correction, not a
    bare SSA friction coefficient. The DIVA and Blatter-Pattyn runs instead start from the
    bare `beta` and derive their own vertical-shear response from Pagos' own 3D viscosity and
    geometry — DIVA via [`diva_update!`](@ref)'s `β_eff`, BP by resolving the column directly
    and reading `beta_eff` (BP has no `F₂` correction to derive, `roadmaps/blatter-pattyn.md`
    §2.2) — which will not exactly reproduce Yelmo's `beta_eff`, not least because of the
    layering mismatch noted in `yelmo-pagos.jl`. So this comparison is "SSA fed Yelmo's
    DIVA-corrected friction" vs. "DIVA/BP deriving their own correction from the same raw
    inputs" — informative about whether Pagos' higher-order paths do something structurally
    different from SSA at all, but not a controlled ablation of the friction law alone.

Same geometry, same friction, same viscosity magnitude — the only difference is which
momentum balance resolves vertical shear, and how. Each call to [`run_solve`](@ref) leaves
at most one `MechanicState` alive at a time; Blatter-Pattyn's is by far the largest (nine
column tensor fields against SSA/DIVA's handful of 2D ones), so this script never holds two
solves' `MechanicState`s at once.
=#
include(joinpath(@__DIR__, "helpers.jl"))

ssa  = run_solve(SSAMomentumBalance(), grid, rt, mask; SOLVER_KWARGS...)
println("SSA:  ", (; ssa.converged, ssa.iterations, ssa.elapsed, ssa.residual))

diva = run_solve(DIVAMomentumBalance(), grid, rt, mask; SOLVER_KWARGS...)
println("DIVA: ", (; diva.converged, diva.iterations, diva.elapsed, diva.residual))

bp   = run_solve(BlatterPattynMomentumBalance(), grid, rt, mask; SOLVER_KWARGS...)
println("BP:   ", (; bp.converged, bp.iterations, bp.elapsed, bp.residual))

#=
## Pairwise diagnostics

`diff(a, b) = a - b`, on-ice only, with the usual RMSE / mean|Δ| / median relative
difference / max|Δ| summary.
=#

function pairwise_stats(name, a, b)
    d = on_ice(a .- b)
    vals = filter(!isnan, d)
    rmse = sqrt(mean(abs2, vals))
    rel  = vals ./ max.(filter(!isnan, on_ice(b)), 1.0)   # avoid /0 on stagnant ice
    @printf("%s (on-ice, %d cells): RMSE = %.4g m/yr, mean|Δ| = %.4g m/yr, median rel. diff = %.4g%%, max|Δ| = %.4g m/yr\n",
           name, length(vals), rmse, mean(abs, vals), 100 * median(abs.(rel)), maximum(abs, vals))
    return d
end

diff_bp_ssa  = pairwise_stats("BP vs SSA",   bp.speed,   ssa.speed)
diff_bp_diva = pairwise_stats("BP vs DIVA",  bp.speed,   diva.speed)
diff_diva_ssa = pairwise_stats("DIVA vs SSA", diva.speed, ssa.speed)

#=
## Figure

Top row (speed, shared colour range): SSA, DIVA, BP. Bottom row (Δ speed, shared symmetric
colour range): BP - SSA, BP - DIVA, DIVA - SSA.
=#

fig1 = Figure(size = (1650, 1120))

crange_top = (0, quantile(filter(!isnan, on_ice(bp.speed)), 0.995))
for (col, (data, title)) in enumerate((
    (on_ice(ssa.speed),  "Pagos SSA"),
    (on_ice(diva.speed), "Pagos DIVA"),
    (on_ice(bp.speed),   "Pagos Blatter-Pattyn"),
))
    ax = Axis(fig1[1, col], xlabel = "x (km)", ylabel = col == 1 ? "y (km)" : "",
        aspect = DataAspect(), title = title)
    hm = heatmap!(ax, xc, yc, data; colorrange = crange_top)
    col == 3 && Colorbar(fig1[1, 4], hm, label = "speed (m/yr)")
end

drange_bottom = maximum(filter(!isnan, abs.(vcat(diff_bp_ssa[:], diff_bp_diva[:], diff_diva_ssa[:]))))
for (col, (data, title)) in enumerate((
    (diff_bp_ssa,   "BP - SSA"),
    (diff_bp_diva,  "BP - DIVA"),
    (diff_diva_ssa, "DIVA - SSA"),
))
    ax = Axis(fig1[2, col], xlabel = "x (km)", ylabel = col == 1 ? "y (km)" : "",
        aspect = DataAspect(), title = title)
    hm = heatmap!(ax, xc, yc, data; colorrange = (-drange_bottom, drange_bottom),
        colormap = cgrad(:RdBu, rev = true))
    col == 3 && Colorbar(fig1[2, 4], hm, label = "Δ speed (m/yr)")
end

save("$figdir/ssa-diva-bp.png", fig1)
fig1

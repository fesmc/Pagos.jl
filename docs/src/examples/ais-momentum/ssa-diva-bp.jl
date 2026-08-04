#=

# SSA vs. DIVA

!!! warning "Not quite apples-to-apples on friction"
    The SSA run uses Yelmo's own `beta_eff` directly but that field
    is *already* Yelmo's own DIVA vertical-shear correction, not a
    bare SSA friction coefficient. The DIVA run instead starts from the bare `beta` and
    derives its own `β_eff` via [`diva_update!`](@ref), from Pagos' own 3D viscosity and
    geometry — which will not exactly reproduce Yelmo's `beta_eff`, not least because of the
    layering mismatch noted in `yelmo-pagos.jl`. So this comparison is "SSA fed Yelmo's
    DIVA-corrected friction" vs. "DIVA deriving its own correction from the same raw inputs" —
    informative about whether Pagos' DIVA path does something structurally different from SSA
    at all, but not a controlled ablation of the friction law alone.

Same geometry, same friction, same viscosity magnitude — the only difference is whether the
momentum balance resolves vertical shear. Each call to [`run_solve`](@ref) leaves at most
one `MechanicState` alive at a time.
=#
include(joinpath(@__DIR__, "helpers.jl"))

ssa  = run_solve(SSAMomentumBalance(), grid, rt, mask; SOLVER_KWARGS...)
println("SSA:  ", (; ssa.converged, ssa.iterations, ssa.elapsed, ssa.residual))

diva = run_solve(DIVAMomentumBalance(), grid, rt, mask; SOLVER_KWARGS...)
println("DIVA: ", (; diva.converged, diva.iterations, diva.elapsed, diva.residual))

on_ice_mask = H_ice .> 0
diff = on_ice(diva.speed .- ssa.speed)
diff_vals = filter(!isnan, diff)
rmse = sqrt(mean(abs2, diff_vals))
rel  = diff_vals ./ max.(filter(!isnan, on_ice(ssa.speed)), 1.0)   # avoid /0 on stagnant ice

@printf("SSA vs DIVA (on-ice, %d cells): RMSE = %.4g m/yr, mean|Δ| = %.4g m/yr, median rel. diff = %.4g%%, max|Δ| = %.4g m/yr\n",
       length(diff_vals), rmse, mean(abs, diff_vals), 100 * median(abs.(rel)), maximum(abs, diff_vals))

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
hm3 = heatmap!(ax3, xc, yc, diff; colorrange = (-drange, drange), colormap = cgrad(:RdBu, rev=true))
Colorbar(fig1[1, 5], hm3, label = "Δ speed (m/yr)")
save("$figdir/ssa-diva.png", fig1)
fig1

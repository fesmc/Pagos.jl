#=

# DIVA: `FixedTuning()` defaults vs. `AutotunedDynamicRelaxation()`

The DIVA counterpart of `roadmaps/PT-autotune.md` Phase 2's headline SSA result (6.3× on
this exact geometry, `gamma = 0.2` found by a hand scan there). No DIVA hand-tuned value
exists anywhere to reuse, and a genuine hand scan means several full solves before the
comparison even starts (per the earlier discussion) — so this compares the **untuned
library default** `FixedTuning()` (`theta_v = 0.6, gamma = 1`) against the autotuner,
rather than simulating a hand search.
=#
include(joinpath(@__DIR__, "helpers.jl"))

diva = run_solve(DIVAMomentumBalance(), grid, rt, mask; SOLVER_KWARGS...)
println("DIVA, AutotunedDynamicRelaxation(): ",
       (; diva.converged, diva.iterations, diva.elapsed, diva.residual))

diva_fixed = run_solve(DIVAMomentumBalance(), grid, rt, mask;
                       (; SOLVER_KWARGS..., tuning = FixedTuning(), maxiter = 1000)...)
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
save("$figdir/hand-auto-tuned.png", fig3)
fig3

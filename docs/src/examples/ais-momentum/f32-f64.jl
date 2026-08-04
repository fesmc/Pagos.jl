#=

# DIVA: Float64 vs. Float32

A genuinely Float32 [`StaggeredGrid`](@ref) — passing `T = Float32` to [`run_solve`](@ref)
alone would *not* do this (see its docstring note); `MechanicState`'s field types come from
the grid, so a separate grid is what actually changes the arithmetic. `_sigma_axis` casts
`layering`'s values to whichever `T` the grid asks for, so the same `Float64`-built
`layering` object is reused rather than needing its own Float32 copy.
=#
include(joinpath(@__DIR__, "helpers.jl"))

diva = run_solve(DIVAMomentumBalance(), grid, rt, mask; SOLVER_KWARGS...)
println("DIVA Float64: ", (; diva.converged, diva.iterations, diva.elapsed, diva.residual))

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
       "by the DIVA solve above, Float32 ones compile here for the first time — some ",
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
save("$figdir/f32-f64.png", fig2)
fig2

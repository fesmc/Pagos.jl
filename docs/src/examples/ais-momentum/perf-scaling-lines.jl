#=

# Performance scaling: wall time vs. resolution

Renders `perf-scaling.csv` (written by `perf-scaling.jl`) as wall time against grid spacing
for the four legs of the sweep. Colour separates the two solvers — linear direct solve
`Cycled(1)`, Gershgorin PT `Cycled(2)` — and linestyle separates the two backends, CPU dashed
and GPU solid, so a leg is identified by two independent visual channels rather than four
arbitrary colours.

Both axes are logarithmic and the x-axis is reversed, so refinement runs left → right and a
power law `t ∝ Δx^(-p)` shows up as a straight line of slope `p`. The grey guide is the ideal
`t ∝ n_dof ∝ Δx^(-2)`: a leg steeper than it is scaling worse than the growth of the problem.

```
julia --project=docs docs/src/examples/ais-momentum/perf-scaling-lines.jl
```
=#
using Pkg
Pkg.activate(joinpath(@__DIR__, "../../.."))
using CairoMakie, DelimitedFiles, Printf

figdir = joinpath(@__DIR__, "figs")
mkpath(figdir)
set_theme!(theme_latexfonts())

#=
## Load the sweep

Rows arrive in invocation order (the 4 km leg runs in its own process and appends), so sort
by resolution rather than assuming the file order. A leg that never ran is `NaN` and is
dropped per-series instead of breaking the whole line.
=#
csv_path = joinpath(@__DIR__, "perf-scaling.csv")
raw, header = readdlm(csv_path, ','; header = true)
col(name) = Float64.(raw[:, findfirst(==(String(name)), vec(header))])

perm = sortperm(col(:resolution_km); rev = true)   # coarse → fine
res  = col(:resolution_km)[perm]
dof  = col(:n_dof)[perm]

series = (
    (:elapsed_lin_cpu, "Linear",  "linear solve",     "CPU", 1, :dash),
    (:elapsed_lin_gpu, "Linear",  "linear solve",     "GPU", 1, :solid),
    (:elapsed_pt_cpu,  "PT",      "psuedo-transient", "CPU", 2, :dash),
    (:elapsed_pt_gpu,  "PT",      "psuedo-transient", "GPU", 2, :solid),
)
marker_of(backend) = backend == "CPU" ? :utriangle : :circle

#=
## The figure

`Cycled(i)` pulls colour `i` from the theme's palette without hard-coding a hex value, which
keeps the two solvers consistent with every other figure drawn under the same theme.

The figure is built from whichever subset of `series` it is handed, so the linear-only variant
below is the same plot with two legs removed rather than a second, drifting copy of the code.
Both legend groups are derived from the legs actually drawn, so dropping the PT legs drops the
PT legend entry with them.
=#
function scaling_figure(legs)
    fig = Figure(size = (700, 600), fontsize = 22)
    ax = Axis(fig[1, 1],
        # title = "DIVA momentum solve: wall time vs. resolution",
        xlabel = "Horizontal resolution (km)",
        ylabel = "Wall time (s)",
        xscale = log2, yscale = log10, xreversed = true,
        xticks = (res, string.(Int.(res))),
        yminorticksvisible = true, yminorgridvisible = true,
        yminorticks = IntervalsBetween(9),
    )

    for (colname, _, _, backend, cyc, linestyle) in legs
        t = col(colname)[perm]
        ok = isfinite.(t)
        any(ok) || continue
        scatterlines!(ax, res[ok], t[ok];
            color = Cycled(cyc), linestyle, linewidth = 3,
            marker = marker_of(backend), markersize = 20)
    end

    # Ideal `t ∝ n_dof`, anchored at the coarsest linear-CPU point so it is read as a slope and
    # not as a competing measurement.
    # t_ref = col(:elapsed_lin_cpu)[perm][1] .* dof ./ dof[1]
    # lines!(ax, res, t_ref; color = (:gray, 0.7), linestyle = :dot, linewidth = 2)

    #=
    Two small legends rather than one four-entry legend: the reader looks up "which solver" and
    "which backend" separately, which is exactly how the encoding is built.
    =#
    solvers  = unique([(leg[5], leg[3]) for leg in legs])   # (colour index, label)
    backends = unique([(leg[4], leg[6]) for leg in legs])   # (name, linestyle)

    solver_entries = [LineElement(color = Cycled(c), linewidth = 3) for (c, _) in solvers]
    backend_entries = [
        [LineElement(color = :black, linewidth = 3, linestyle = ls),
         MarkerElement(color = :black, marker = marker_of(b), markersize = 20)]
        for (b, ls) in backends
    ]
    ref_entries = [LineElement(color = (:gray, 0.7), linewidth = 2, linestyle = :dot)]

    axislegend(ax,
        [solver_entries, backend_entries], # ref_entries
        [[lab for (_, lab) in solvers], [b for (b, _) in backends]], # [L"$t \propto n_\mathrm{dof}$"]],
        ["solver", "backend"]; # , "reference"
        position = :lt, labelsize = 14, titlesize = 14, patchsize = (30, 14),
        # framevisible = false,
        rowgap = 2, titlegap = 2, groupgap = 8,
    )
    ylims!(ax, 1e-2, 2e2)
    return fig
end

fig = scaling_figure(series)
ylims!(ax, 1e-2, 2e2)
save(joinpath(figdir, "perf-scaling-lines.png"), fig)

#=
## Linear solve alone

The PT legs span two orders of magnitude more than the linear ones and set the y-range, which
squashes the CPU/GPU crossover that the direct solve is actually being read for. Same figure,
same colour, PT dropped — so the linear pair gets the full height of the axis.
=#
fig_lin = scaling_figure(filter(leg -> leg[2] == "Linear", collect(series)))
# ax = axes(fig_lin[1, 1])
# ylims!(ax, 1e-2, 2e2)
save(joinpath(figdir, "perf-scaling-lines-linear.png"), fig_lin)

#=
## Effective scaling exponents

`p` in `t ∝ Δx^(-p)` from a least-squares fit through all finite points of a leg; `p = 2` is
the ideal (cost proportional to the number of degrees of freedom).
=#
for (colname, solver, _, backend, _, _) in series
    t = col(colname)[perm]
    ok = isfinite.(t)
    count(ok) ≥ 2 || continue
    x = log.(res[ok])
    y = log.(t[ok])
    p = -((x .- sum(x) / length(x))' * (y .- sum(y) / length(y))) / sum(abs2, x .- sum(x) / length(x))
    @printf("%-7s %-3s  p = %.2f   (%.3f s → %.3f s)\n", solver, backend, p, t[ok][1], t[ok][end])
end

fig

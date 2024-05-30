using NCDatasets, CairoMakie

function load_runtimes(file::String)
    ds = NCDataset(file, "r")
    dx = ds["resolution"][:]
    runtime = ds["runtime"][:]
    close(ds)
    return dx, runtime
end

bmdir = "$(@__DIR__)/../.."
dx_pt, runtime_pt = load_runtimes("$bmdir/data/PTloop-runtimes.nc")
dx_lin, runtime_lin = load_runtimes("$bmdir/data/linearsolve-runtimes.nc")
runtime_lin[end] = NaN
powervec = log2.(dx_pt ./ 1e3)

fig = Figure(size = (800, 600), fontsize = 20)
ax = Axis(fig[1, 1], yscale = log10, yminorticks = IntervalsBetween(9))
ax.xlabel = "Resolution (km)"
ax.ylabel = "Runtime (s)"
ax.title = "Pseudo-transient vs. linear solve"

lines!(ax, log2.(dx_pt ./ 1e3), runtime_pt, linewidth = 2, label = "Pseudo-transient")
lines!(ax, log2.(dx_lin ./ 1e3), runtime_lin, linewidth = 2, label = "Linear solve")
axislegend(ax, position = :rt)
ax.xticks = (powervec, ["$dx" for dx in Int.(2 .^ powervec)])
ax.yminorticksvisible = true
ax.yminorgridvisible = true
fig

save("$bmdir/plots/log-linear-vs-PT.png", fig)
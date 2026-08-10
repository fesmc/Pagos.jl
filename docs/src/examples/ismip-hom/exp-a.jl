#=

# ISMIP-HOM A: Pagos' Blatter-Pattyn against a full-Stokes ensemble

`pagos-roadmaps/blatter-pattyn.md` Phase 3. Experiment **A** — the 3D one: a sinusoidal bed
bumpy in both `x` and `y`, no slip, periodic, at domain lengths
`L = 160, 80, 40, 20, 10 km`. Experiment B lives in `exp-b.jl`; everything shared between
the two is in `helpers.jl`, including the geometry, the file format, the `Δp` sign
convention warning, and why `L = 5 km` is excluded.

A is the harder of the pair and the one that cannot be faked: the transverse bed variation
is exactly what a flowline model has no representation for, so it is where a first-order
solver's handling of the *membrane* terms — the `∂/∂y` half of the stress divergence, which
B never exercises — is actually under test.
=#
include(joinpath(@__DIR__, "helpers.jl"))

RUN_PAGOS = @isdefined(RUN_PAGOS) ? RUN_PAGOS : true

LABELS = model_label.(MODELS)
A = [[load_expA(d, L) for L in LENGTHS_KM] for d in MODELS]

println("full-Stokes ensemble, experiment A:")
for (m, lab) in enumerate(LABELS)
    a = A[m][1]
    @printf("  %-5s — %d×%d, %s-fastest, slice at y/L = %.4f (nearest %.2f)\n",
            lab, a.n, a.n, a.x_fastest ? "x" : "y", a.y_slice, SLICE_Y)
end

## Grouped by domain length rather than by model — that is the shape `fs_envelope` consumes,
## since the band is taken across the ensemble at one `L`.
refsA = [[A[m][i].slice for m in eachindex(MODELS)] for i in eachindex(LENGTHS_KM)]

#=
## The `Δp` sign convention across the ensemble

Detected, not assumed — see the warning in `helpers.jl`. With five submissions this stops
being a curiosity: one member with the opposite convention would double the apparent width
of the `Δp` band and put its centre near zero, which would read as "full Stokes cannot agree
on the pressure anomaly" rather than "one file writes `p_H − p_I`".
=#
flipsA = [[dp_orientation(p) for p in refsA[i]] for i in eachindex(LENGTHS_KM)]
dpmaskA = dp_ensemble_mask(refsA)
println()
for (m, lab) in enumerate(LABELS)
    @printf("  %-5s — Δp: %-38s %s\n", lab,
            flipsA[end][m] > 0 ? "compressive-positive on the stoss side" :
                                 "OPPOSITE SIGN (normalised)",
            dpmaskA[m] ? "" : "— EXCLUDED (fails the shallow-ice limit)")
end

#=
## Running Pagos

`USE_GPU` picks the backend; `helpers.jl` explains why the state is filled on the host and
adapted in bulk. Set `RUN_PAGOS = false` before `include`ing to plot the reference ensemble
alone.
=#
pagosA = nothing
if RUN_PAGOS
    println("\nPagos Blatter-Pattyn (periodic halo, hand-damped PT) on ",
            USE_GPU ? "GPU" : "CPU", ":")
    pagosA = Vector{Any}(undef, length(LENGTHS_KM))
    for (i, L) in enumerate(LENGTHS_KM)
        r = run_pagos_bp(:A, L); G = pagos_fields(r)
        j = argmin(abs.(G.y .- SLICE_Y))
        pagosA[i] = (; n = r.nx, G.x, G.y, G.vx, G.vy, G.vz, G.txz, G.tyz, G.dp,
                       slice = (; G.x, vx = G.vx[:, j], vz = G.vz[:, j],
                                  txz = G.txz[:, j], dp = G.dp[:, j]),
                       y_slice = G.y[j], x_fastest = true, r.converged, r.err)
        @printf("  L = %3d km   nx = %d, nz = %d   %7.1fs  %6d iters  err = %.1e %s\n",
                L, r.nx, r.nz, r.elapsed, r.iters_used, r.err, r.converged ? "ok" : "!!")
    end
end
pagos_slices = RUN_PAGOS ? [p.slice for p in pagosA] : nothing

#=
`helpers.jl` explains why a hand-rolled periodic loop is needed here too, and why it can
reuse the library's own [`pseudo_rate!`](@ref) unmodified where BP's could not: DIVA's
per-iteration step has no terrain-following correction to interleave a halo refresh into.
Solved with the default [`GershgorinPseudoTimeStep`](@ref) tuning of `Δτ`, the same bound
`run_pagos_bp` uses, extended with the basal-drag term `diva_update!`'s `β_eff` feeds it.
=#
pagosA_diva = nothing
if RUN_PAGOS
    println("\nPagos DIVA (periodic halo, Gershgorin-tuned PT) on ",
            USE_GPU ? "GPU" : "CPU", ":")
    pagosA_diva = Vector{Any}(undef, length(LENGTHS_KM))
    for (i, L) in enumerate(LENGTHS_KM)
        r = run_pagos_diva(:A, L); G = pagos_fields(r)
        j = argmin(abs.(G.y .- SLICE_Y))
        pagosA_diva[i] = (; n = r.nx, G.x, G.y, G.vx, G.vy, G.vz, G.txz, G.tyz, G.dp,
                            slice = (; G.x, vx = G.vx[:, j], vz = G.vz[:, j],
                                       txz = G.txz[:, j], dp = G.dp[:, j]),
                            y_slice = G.y[j], x_fastest = true, r.converged, r.err)
        @printf("  L = %3d km   nx = %d, nz = %d   %7.1fs  %6d iters  err = %.1e %s\n",
                L, r.nx, r.nz, r.elapsed, r.iters_used, r.err, r.converged ? "ok" : "!!")
    end
end
pagos_diva_slices = RUN_PAGOS ? [p.slice for p in pagosA_diva] : nothing

#=
## Profiles at `ŷ ≈ 0.25`

The shaded band is the full-Stokes min–max envelope at each domain length; the solid line is
Pagos BP, the dashed line Pagos DIVA. The slice is where `sin(ω y) = 1`, so the bed bumps are
cut at their full 500 m amplitude.
=#
figA = profile_figure(refsA, pagos_slices,
                      "ISMIP-HOM A — bumpy bed, no slip; slice at y/L ≈ 0.25",
                      length(MODELS); flips = flipsA, dpmask = dpmaskA,
                      pagos2 = pagos_diva_slices)
save(joinpath(figdir, "ismip-hom-a.png"), figA)
figA

#=
## How the solution scales with `L`
=#
figAscale = scaling_figure(refsA, pagos_slices,
                           "ISMIP-HOM A — dependence on domain length"; dpmask = dpmaskA)
save(joinpath(figdir, "ismip-hom-a-scaling.png"), figAscale)
figAscale

#=
## Experiment A is genuinely 3D — the slice hides that

A single `ŷ` slice cannot show that A's bed varies in `y` as well, which is the whole
difference from B. Surface `vx` over the full normalized domain, one panel per `L`, each with
its own colour range because the magnitudes span roughly 7× between `L = 160 km` and
`L = 10 km`. The dashed line marks the profile slice above.

One figure per model rather than a difference map: the submissions use different grids (41×41
cell centres against 61×61 nodes), and interpolating one onto the other to subtract would
manufacture a field neither model produced. The band in the profile figure is where the
ensemble is compared quantitatively, and it says so.
=#
function maps_figure(series, label)
    fig = Figure(size = (1500, 460))
    Label(fig[0, 1:5], "ISMIP-HOM A ($label) — surface velocity vx(zs), m/yr";
          fontsize = 19, font = :bold, padding = (0, 0, 6, 0))
    for (i, (L, a)) in enumerate(zip(LENGTHS_KM, series))
        gl = GridLayout(fig[1, i])
        ax = Axis(gl[1, 1]; title = "L = $L km", xlabel = "normalized x",
                  ylabel = i == 1 ? "normalized y" : "")
        hm = heatmap!(ax, a.x, a.y, a.vx; colormap = :viridis)
        hlines!(ax, [a.y_slice]; color = (:white, 0.85), linestyle = :dash, linewidth = 1.5)
        Colorbar(gl[1, 2], hm; width = 10, ticklabelsize = 10)
        ## `Aspect` on the *cell* rather than `aspect = DataAspect()` on the axis: the latter
        ## shrinks the axis inside a cell that keeps its full width, stranding the colorbar a
        ## panel-width away. This constrains the cell itself to square, so the two sit adjacent.
        colsize!(gl, 1, Aspect(1, 1.0))
        colgap!(gl, 8)
    end
    return fig
end

for (m, lab) in enumerate(LABELS)
    save(joinpath(figdir, "ismip-hom-a-maps-$lab.png"), maps_figure(A[m], lab))
end
RUN_PAGOS && save(joinpath(figdir, "ismip-hom-a-maps-pagos.png"),
                  maps_figure(pagosA, "Pagos BP"))
maps_figure(A[1], LABELS[1])

#=
## Numbers behind the figures
=#
results_table(:A, refsA, pagos_slices, LABELS; dpmask = dpmaskA)

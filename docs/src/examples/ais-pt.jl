#=

# [Antarctic pseudo-transient solve from real geometry](@id ais_pt)

This example exercises the Chmy-native, C-grid staggered pseudo-transient (PT) momentum
solver ([`pseudo_transient!`](@ref), Sandip et al. 2024) on real Antarctic Ice Sheet (AIS)
geometry, rather than an analytic or idealized test case. It is a code-state check, not a
validated simulation: the goal is to see what the current solver produces on a realistic
761×761, 8 km domain, and it doubles as the real-geometry measurement of the Duretz et al.
(2026) autotuner (`roadmaps/PT-autotune.md`, Phase 2) — see the solver settings below.

Fields are read from a Yelmo restart file: `x`/`y` coordinates, bed/surface topography and
ice thickness, the depth-averaged viscosity `visc_bar` and the effective basal friction
coefficient `beta_eff`. Viscosity is a *fixed* input (no Glen-law viscosity continuation);
the basal stress is **not** — it is recomputed as `β_eff · u` every PT iteration, which is
what makes the momentum balance a well-posed problem for `u` (see the warning below).

Velocity starts at zero everywhere and is relaxed in pseudo-time until either the
nondimensional momentum residual converges or a maximum iteration count is reached. The
converged field is plotted next to Yelmo's own `uxy_bar` from the same restart, which is the
reference this example is checked against — the aim is the right pattern and magnitude, not
a point-by-point match (Yelmo solves DIVA with its own discretization, boundary and
calving-front treatment; Pagos is in the SSA limit with a placeholder `Neumann(0)` ring).

!!! warning "Why friction is solved for, not prescribed"
    An earlier version of this example wrote Yelmo's `taub_acx`/`taub_acy` straight into
    `mech.stress.base_x`/`base_y` and held them fixed ([`NoFrictionUpdate`](@ref)),
    bypassing the friction law. That is not merely an approximation, it removes the only
    velocity-dependent restoring force in the balance: with `τ_b` a constant field the
    operator has no zeroth-order term at all, the problem inherits a rigid-translation null
    space that the `Neumann(0)` ring does not pin, and the iteration drifts linearly forever
    instead of converging (the residual plateaus while `max|u|` grows without bound). Using
    `beta_eff` costs nothing — it is in the same restart file — and makes the problem the
    solver is actually built for.

!!! warning "Known simplifications"
    - **Face staggering.** Yelmo stores `acx`/`acy` fields with the same shape as `aa`
      fields (`taub_acx[i,j]` is the face at `i+1/2`); Pagos' `StaggeredGrid` uses an
      explicit `Vertex` axis with one more point than `Center` (vertex `j` sits between
      centres `j-1` and `j`, see `roadmaps/yelmo-staggering.md`). Here every input is a
      cell-centred (`aa`) field, so no such shift is needed — `beta_eff` is staggered onto
      the velocity faces by Pagos itself, inside [`basalstress!`](@ref).
    - **Units.** `Constants` is SI throughout, so `visc_bar` (Yelmo units: Pa yr) and
      `beta_eff` (Pa yr m⁻¹) are converted to Pa s and Pa s m⁻¹; the solve runs in SI (m/s)
      and the converged velocity is converted back to m/yr only for reporting/plotting.
    - **SSA, not DIVA.** `beta_eff` carries Yelmo's DIVA vertical-shear correction, but
      Pagos applies it in the SSA limit (`F₂` is Phase 3 work), so the two do not solve
      quite the same equation.

=#

using Pagos, NCDatasets, CairoMakie, Statistics

restart_file = "/home/jan/pCloudSync/PhD/Projects/Ice-Sheet-Modelling/ice-data-pagos/yelmo_restart_ais_8km.nc"

#=
## Load the restart fields

Everything below is squeezed to 2D (the file carries a size-1 `time` dimension) and cast
to `Float64`, matching the grid's element type. `uxy_bar` is Yelmo's own converged
depth-averaged speed — the reference, never an input to the solve.
=#

f64(x) = Float64.(dropdims(x; dims = ndims(x)))

xc, yc, H_ice, z_srf, z_bed, visc_bar, beta_eff, f_grnd, uxy_bar =
NCDataset(restart_file) do ds
    (Float64.(ds["xc"][:]), Float64.(ds["yc"][:]),
     f64(ds["H_ice"][:, :, :]), f64(ds["z_srf"][:, :, :]), f64(ds["z_bed"][:, :, :]),
     f64(ds["visc_bar"][:, :, :]), f64(ds["beta_eff"][:, :, :]),
     f64(ds["f_grnd"][:, :, :]), f64(ds["uxy_bar"][:, :, :]))
end

nx, ny = length(xc), length(yc)
dx = (xc[2] - xc[1]) * 1e3   # km -> m
dy = (yc[2] - yc[1]) * 1e3
lx, ly = nx * dx, ny * dy

#=
## Build the grid, mechanics state and ice mask

`StaggeredGrid`'s cell-centre convention (`-lx/2 + dx/2 : dx : lx/2 - dx/2`) reproduces the
restart file's `xc`/`yc` exactly for this `(nx, dx)`, so `aa`-node fields can be copied over
index-for-index with no reprojection.

The mask handed to the solver is [`momentum_mask!`](@ref)'s `is_momentum_solved`, **not**
`is_ice`. Antarctica's ice at this resolution includes 48 cells (6 patches, 0.02% of the
ice) that are detached from the grounded sheet. An iceberg has no basal drag and no membrane
connection, so nothing balances its driving stress: the momentum balance restricted to it is
singular in its rigid-translation modes, and no solver setting can converge a problem with
no answer. Under a max-norm stopping criterion those few faces set `err` for the other
211258 cells — with `is_ice` this solve runs to `maxiter` with the residual pinned at 3.8e-3
and a berg at ~6000 m/yr; with `is_momentum_solved` it converges. Advection would still be
given `IceMask(is_ice, is_ice_neighbour)`, so bergs keep advecting and calving; they are
only excused from a force balance they cannot satisfy.
=#

grid = StaggeredGrid(Float64, lx, ly, dx, dy)
rt   = Runtime(grid)
mech = MechanicState(grid)
topo = TopographicState(grid)
cst  = Constants{Float64}()

#=
## Load geometry and material fields

`topography.thickness`/`surface`, `material.viscosity_depthaveraged` and
`friction.beta_eff` are read by kernels that touch their halo (the membrane-stress corner
viscosity in particular), so they are filled with `fill_from_grid!`, a discrete analogue of
the test helpers' `fill_analytic!`: the domain's own edge values are extrapolated
(zero-gradient) into the halo ring rather than left at the allocation default.

Off-ice viscosity is set to the domain's maximum rather than left at the file's `0`. With
[`GershgorinPseudoTimeStep`](@ref) below this value is never actually read — every
coefficient of that bound is evaluated through `mask`, so an ice-free cell contributes the
same zero it contributes to the residual — but it keeps the state well-defined for the
[`ViscosityPseudoTimeStep`](@ref) path, which divides by an interpolated viscosity
unconditionally and would produce `Inf` on a true `0`.
=#

function fill_from_grid!(f, data)
    ni, nj = size(data)
    for k in axes(interior(f), 3), j in -1:(nj + 2), i in -1:(ni + 2)
        ic, jc = clamp(i, 1, ni), clamp(j, 1, nj)
        f[i, j, k] = data[ic, jc]
    end
    return f
end

spy            = cst.seconds_per_year
visc_bar_si    = visc_bar .* spy                    # Pa yr     -> Pa s
beta_eff_si    = beta_eff .* spy                    # Pa yr m⁻¹ -> Pa s m⁻¹
visc_off_ice   = maximum(visc_bar_si)
visc_for_solve = ifelse.(H_ice .> 0, visc_bar_si, visc_off_ice)

fill_from_grid!(mech.topography.thickness, H_ice)
fill_from_grid!(mech.topography.surface, z_srf)
fill_from_grid!(mech.material.viscosity_depthaveraged, visc_for_solve)
fill_from_grid!(mech.friction.beta_eff, beta_eff_si)

# The mask: `icemasks!` derives `is_ice` from the thickness, `momentum_mask!` grows
# `is_momentum_solved` out from the grounded seed through connected ice. `is_grounded` is
# read, not derived — here from Yelmo's sub-grid grounded fraction.
fill_from_grid!(topo.thickness.ice, H_ice)
fill_from_grid!(topo.mask.is_grounded, f_grnd .> 0)
icemasks!(topo, rt)
momentum_mask!(topo, rt)
mask = IceMask(topo.mask.is_momentum_solved)

n_detached = count(asarray(topo.mask.is_ice) .& .!asarray(topo.mask.is_momentum_solved))
@show n_detached

# Velocity starts at zero (the `MechanicState` allocation default; set explicitly here for
# clarity since this is the whole point of the exercise).
setdata!(mech.velocity.x, 0.0)
setdata!(mech.velocity.y, 0.0)

#=
## Run the solve

Three solver choices matter here. The first is now the default *because* of this example;
the other two still have to be asked for:

 1. **[`GershgorinPseudoTimeStep`](@ref)** rather than [`ViscosityPseudoTimeStep`](@ref).
    Sandip's `Δτ ∝ 1/η` bounds only the membrane part of the operator. It omits the basal
    drag `β/(ρH)`, which is the dominant eigenvalue under grounded Antarctic ice
    (`β_eff` reaches ~5·10¹³ Pa s m⁻¹ here) — with friction solved for rather than
    prescribed, that `Δτ` diverges within ~20 iterations. The Gershgorin bound is built
    from the coefficients the residual kernels actually use, drag included, and is
    mask-aware at the margin so an ice-free neighbour no longer throttles `Δτ` on the very
    faces that carry the calving-front forcing. Only `cfl` is set explicitly below (`0.99`
    rather than the default `0.9`, since the autotuner wants a tight one).
 2. **[`ScaledResidual`](@ref)** instead of [`VelocityIncrement`](@ref). An ice shelf has no
    basal drag, so it relaxes diffusively and its per-iteration velocity increment is
    orders of magnitude below the grounded ice's *from the first iteration*. Any `abstol`
    the grounded ice has to work for is one the shelves satisfy at ~0 velocity, so the
    increment criterion returns `converged = true` with Ross and Ronne empty. `abstol` is
    dimensionless here: the fraction of the driving-stress forcing left unbalanced.
 3. **[`AutotunedDynamicRelaxation`](@ref)** instead of a hand-set `gamma`. This is the one
    place where the paper claim is directly measurable on real geometry. An earlier version
    of this example carried `theta_v = 1.0, gamma = 0.2`, a value found by scanning, and
    converged in 1500 iterations (~140 s). Deriving the damping instead — Gershgorin
    `λ_max`, Rayleigh-quotient `λ_min`, re-estimated every 20 iterations — converges the
    *same* problem to the *same* tolerance in **240 iterations (~24 s), a 6.3× speedup**,
    and the derived `γ ≈ 0.032` shows the hand scan had been ~6× too large. Nothing about
    the scan was careless; the point is that the right value is a property of this mesh and
    this viscosity field, and reading it off the operator beats guessing it.

The fourth ingredient is the iceberg mask built above — without it this same solver runs out
of iterations with the residual pinned by six detached patches. Note that no `theta_v` or
`gamma` appears anywhere: those live on [`FixedTuning`](@ref), the tuning this solve is not
using.
=#

solver = PseudoTransientSolver(grid;
    abstol = 1e-3,          # dimensionless (ScaledResidual)
    maxiter = 2000,
    ncheck = 20,
    printout_every = 50,
    pseudo_timestep = GershgorinPseudoTimeStep(cfl = 0.99),
    convergence = ScaledResidual(),
    friction_update = ActiveFrictionUpdate(),
    tuning = AutotunedDynamicRelaxation(),
)
momentum = SSAMomentumBalance()

result = pseudo_transient!(mech, cst, solver, rt, momentum, mask)
@show result

#=
## Visualize the converged velocity field, against Yelmo's own
=#

ux_c = @views (interior(mech.velocity.x)[1:(end - 1), :, 1] .+
               interior(mech.velocity.x)[2:end, :, 1]) ./ 2 .* spy
uy_c = @views (interior(mech.velocity.y)[:, 1:(end - 1), 1] .+
               interior(mech.velocity.y)[:, 2:end, 1]) ./ 2 .* spy
speed = sqrt.(ux_c .^ 2 .+ uy_c .^ 2)

on_ice(a) = ifelse.(H_ice .> 0, a, NaN)
speed_plot = on_ice(speed)
yelmo_plot = on_ice(uxy_bar)
crange = (0, quantile(filter(!isnan, yelmo_plot), 0.995))

set_theme!(theme_latexfonts())
fig = Figure(size = (1150, 620))
for (col, (data, title)) in enumerate(((speed_plot, "Pagos PT (SSA, prescribed β and η)"),
                                       (yelmo_plot, "Yelmo restart (uxy_bar)")))
    ax = Axis(fig[1, col], xlabel = "x (km)", ylabel = col == 1 ? "y (km)" : "",
        aspect = DataAspect(), title = title)
    hm = heatmap!(ax, xc, yc, data; colormap = :viridis, colorrange = crange)
    contour!(ax, xc, yc, H_ice; levels = [0.5], color = :black, linewidth = 0.75)
    col == 2 && Colorbar(fig[1, 3], hm, label = "Speed (m yr⁻¹)")
end
Label(fig[2, 1:3],
    "iterations = $(result.iterations), converged = $(result.converged), " *
    "scaled residual = $(round(result.error, sigdigits = 3)), " *
    "autotuned γ = $(round(result.damping, sigdigits = 3)) " *
    "(λ_min = $(round(result.lambda_min, sigdigits = 3)))",
    fontsize = 12)

save(joinpath(@__DIR__, "ais-pt-velocity.png"), fig)
fig

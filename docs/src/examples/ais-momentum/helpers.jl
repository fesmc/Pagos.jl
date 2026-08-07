#=
Shared setup for the [Antarctic momentum-balance comparisons](@id ais_momentum): data
loading, grid/topography/mask construction, and the common [`run_solve`](@ref) helper. Each
comparison script `include`s this file first, then runs its own solves.

None of the comparisons is a validated simulation: this is a code-state check on real
geometry, not a benchmark paper.

`load2d`/`load3d` read one variable, drop the trailing size-1 `time` dimension, and cast
directly to `T` in a single pass — no Float64 intermediate that then gets thrown away.
=#
using Pkg
Pkg.activate(joinpath(@__DIR__, "../../.."))
using Pagos, NCDatasets, CairoMakie, Statistics, Printf

#=
`resolution_km` picks which Yelmo restart to load — `8` or `16`, matching the two restarts
that actually exist on disk (`ice-data-pagos/yelmo_restart_ais_{8,16}km.nc`). A script that
wants a specific resolution sets it before `include`ing this file; scripts that don't care
get the `8` km default below.
=#
resolution_km = @isdefined(resolution_km) ? resolution_km : 8
restart_file = if resolution_km == 8
    "/home/jan/pCloudSync/PhD/Projects/Ice-Sheet-Modelling/ice-data-pagos/yelmo_restart_ais_8km.nc"
elseif resolution_km == 16
    "/home/jan/pCloudSync/PhD/Projects/Ice-Sheet-Modelling/ice-data-pagos/yelmo_restart_ais_16km.nc"
else
    error("resolution_km must be 8 or 16, got $resolution_km")
end
figdir = joinpath(@__DIR__, "figs")

load2d(ds, name, T) = T.(dropdims(ds[name][:, :, :]; dims = 3))
load3d(ds, name, T) = T.(dropdims(ds[name][:, :, :, :]; dims = 4))

xc, yc, H_ice, z_srf, z_bed, f_grnd, visc_bar, beta, beta_eff, visc3d, zeta =
NCDataset(restart_file) do ds
    (Float64.(ds["xc"][:]), Float64.(ds["yc"][:]),
     load2d(ds, "H_ice", Float16), load2d(ds, "z_srf", Float16), load2d(ds, "z_bed", Float16),
     load2d(ds, "f_grnd", Float16),
     load2d(ds, "visc_bar", Float32), load2d(ds, "beta", Float32), load2d(ds, "beta_eff", Float32),
     load3d(ds, "visc", Float32),
     Float64.(ds["zeta"][:]))
end

nx, ny = length(xc), length(yc)
nz = length(zeta)
dx = (xc[2] - xc[1]) * 1e3
dy = (yc[2] - yc[1]) * 1e3
lx, ly = nx * dx, ny * dy

#=
## Build the grid, ice masks, and momentum mask

One [`StaggeredGrid`](@ref) with a real column (`nz = 11`, [`QuadraticSigmaTransform`](@ref)),
shared by every comparison below: the [`TopographyMasks`](@ref) are built once, since nothing
about them depends on which momentum balance or element type a given solve uses. Only the masks
are needed here — not a full [`TopographicState`](@ref), which would also carry the unused
mass-balance and elevation fields `run_solve` never touches. The `is_momentum_solved` mask
excludes the ~48 detached iceberg cells a force balance cannot be posed on.
=#
T = Float64
layering = CorrectedVerticalLayering(T, QuadraticSigmaTransform(T, nz))
grid = StaggeredGrid(T, lx, ly, dx, dy, layering)
rt   = Runtime(grid)
masks = TopographyMasks(grid)
cst  = Constants{T}()

function fill_from_grid!(f, data)
    ni, nj = size(data)
    for k in axes(interior(f), 3), j in -1:(nj + 2), i in -1:(ni + 2)
        ic, jc = clamp(i, 1, ni), clamp(j, 1, nj)
        f[i, j, k] = data[ic, jc]
    end
    return f
end

function fill_from_grid3d!(f, data)
    ni, nj, nk = size(data)
    for k in 1:nk, j in -1:(nj + 2), i in -1:(ni + 2)
        ic, jc = clamp(i, 1, ni), clamp(j, 1, nj)
        f[i, j, k] = data[ic, jc, k]
    end
    return f
end

thickness_ice = Field(grid.arch, grid.grid2d, (Center(), Center(), Center()), T; halo = 1)
fill_from_grid!(thickness_ice, H_ice)
fill_from_grid!(masks.is_grounded, f_grnd .> 0)
icemasks!(masks, thickness_ice, rt)
momentum_mask!(masks.is_momentum_solved, masks.is_ice, masks.is_grounded, rt)
mask = IceMask(masks.is_momentum_solved)

n_detached = count(asarray(masks.is_ice) .& .!asarray(masks.is_momentum_solved))
@show n_detached

#=
## Units

No conversion needed. Pagos computes in `(m, yr, Pa)` — see [`Constants`](@ref) — and
Yelmo's restart fields are already `Pa yr` / `Pa yr m⁻¹`, so viscosity, friction and
velocity all line up as read.
=#

visc_off_ice   = maximum(visc_bar)
visc3d_off_ice = maximum(visc3d)
@. visc_bar = ifelse(H_ice > 0, visc_bar, visc_off_ice)
@. visc3d   = ifelse(H_ice > 0, visc3d, visc3d_off_ice)

on_ice(a) = ifelse.(H_ice .> 0, a, NaN);

#=
## Run a solve and extract the bare minimum from the solution

Only a small 2D speed field (`nx × ny`, ~2 MB at Float32, not the ~2 GB `mech` it came from)
and a few scalars survive.
=#

#=
`MomentumBalance3D` (Blatter-Pattyn) iterates the full column `velocity.x`/`y`, not
`velocity.depthaverage_x`/`y` — depth-averaging it down to the same `(nx, ny)` speed field
the SSA/DIVA runs produce is a diagnostic reduction the solve itself has no reason to do, so
it lives here rather than in the library. Weighted by the sigma-axis layer thickness `Δζ_k`
— the same weights `depthaverage!` uses for `µ̄` — read at a single representative column
since `CorrectedVerticalLayering`'s `ζ` layout does not depend on `(i, j)`.
=#
function column_depthaverage_speed(mech, rt)
    (; x, y) = mech.velocity
    nz = size(rt.grid, Center())[3]
    wts = reshape([Δz(rt.grid, Center(), 1, 1, k) for k in 1:nz], 1, 1, nz)

    ubar_x = dropdims(sum(interior(x) .* wts; dims = 3); dims = 3)
    ubar_y = dropdims(sum(interior(y) .* wts; dims = 3); dims = 3)

    return Float32.(sqrt.(
        (@views (ubar_x[1:(end - 1), :] .+ ubar_x[2:end, :]) ./ 2) .^ 2 .+
        (@views (ubar_y[:, 1:(end - 1)] .+ ubar_y[:, 2:end]) ./ 2) .^ 2))
end

function run_solve(momentum, grid, rt, mask; solver_kwargs...)
    T = eltype(grid.grid)
    mech = MechanicState(grid)
    cst_T = Constants{T}()

    fill_from_grid!(mech.topography.thickness, T.(H_ice))
    fill_from_grid!(mech.topography.surface, T.(z_srf))
    fill_from_grid!(mech.material.viscosity_depthaveraged, T.(visc_bar))
    fill_from_grid!(mech.friction.beta_eff, T.(beta_eff))
    fill_from_grid!(mech.friction.beta, T.(beta))
    fill_from_grid3d!(mech.material.viscosity, T.(visc3d))

    if momentum isa MomentumBalance3D
        setdata!(mech.velocity.x, zero(T))
        setdata!(mech.velocity.y, zero(T))
        solver = PseudoTransientSolver(grid, momentum; solver_kwargs...)
    else
        setdata!(mech.velocity.depthaverage_x, zero(T))
        setdata!(mech.velocity.depthaverage_y, zero(T))
        solver = PseudoTransientSolver(grid; solver_kwargs...)
        momentum isa DIVAMomentumBalance && diva_update!(mech, solver, rt, mask)
    end

    elapsed = @elapsed result = pseudo_transient!(mech, cst_T, solver, rt, momentum, mask)

    speed = if momentum isa MomentumBalance3D
        column_depthaverage_speed(mech, rt)
    else
        Float32.(sqrt.(
            (@views (interior(mech.velocity.depthaverage_x)[1:(end - 1), :, 1] .+
                     interior(mech.velocity.depthaverage_x)[2:end, :, 1]) ./ 2) .^ 2 .+
            (@views (interior(mech.velocity.depthaverage_y)[:, 1:(end - 1), 1] .+
                     interior(mech.velocity.depthaverage_y)[:, 2:end, 1]) ./ 2) .^ 2))
    end

    mech = nothing
    GC.gc()
    return (; speed, elapsed, iterations = result.iterations, converged = result.converged,
           residual = result.residual)
end

const SOLVER_KWARGS = (
    abstol = 1e-3, maxiter = 2000, ncheck = 20, printout_every = 50,
    pseudo_timestep = GershgorinPseudoTimeStep(cfl = 0.8),
    convergence = ScaledResidual(),
    friction_update = ActiveFrictionUpdate(),
    tuning = AutotunedDynamicRelaxation(cadence = 50),
)

set_theme!(theme_latexfonts())

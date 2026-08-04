#=
Shared setup for the [Antarctic momentum-balance comparisons](@id ais_momentum): data
loading, grid/topography/mask construction, and the common [`run_solve`](@ref) helper. Each
comparison script `include`s this file first, then runs its own solves.

None of the comparisons is a validated simulation — see `ais-pt.jl`'s own framing, which
still applies: this is a code-state check on real geometry, not a benchmark paper.

`load2d`/`load3d` read one variable, drop the trailing size-1 `time` dimension, and cast
directly to `T` in a single pass — no Float64 intermediate that then gets thrown away.
=#
using Pkg
Pkg.activate(joinpath(@__DIR__, "../../.."))
using Pagos, NCDatasets, CairoMakie, Statistics, Printf

restart_file = "/home/jan/pCloudSync/PhD/Projects/Ice-Sheet-Modelling/ice-data-pagos/yelmo_restart_ais_8km.nc"
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
## Build the grid, topography state, and ice mask

One [`StaggeredGrid`](@ref) with a real column (`nz = 11`, [`QuadraticSigmaTransform`](@ref)),
shared by every comparison below: [`TopographicState`](@ref) and the mask are built once, since
nothing about them depends on which momentum balance or element type a given solve uses. The
`is_momentum_solved` mask excludes the ~48 detached iceberg cells a force balance cannot be posed on.
=#

layering = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, nz))
grid = StaggeredGrid(Float64, lx, ly, dx, dy, layering)
rt   = Runtime(grid)
topo = TopographicState(grid)
cst  = Constants{Float64}()

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

fill_from_grid!(topo.thickness.ice, H_ice)
fill_from_grid!(topo.mask.is_grounded, f_grnd .> 0)
icemasks!(topo, rt)
momentum_mask!(topo, rt)
mask = IceMask(topo.mask.is_momentum_solved)

n_detached = count(asarray(topo.mask.is_ice) .& .!asarray(topo.mask.is_momentum_solved))
@show n_detached

#=
## Units

No conversion needed. Pagos computes in `(m, yr, Pa)` — see [`Constants`](@ref) — and
Yelmo's restart fields are already `Pa yr` / `Pa yr m⁻¹`, so viscosity, friction and
velocity all line up as read. This block used to multiply through by `seconds_per_year`
to reach an SI-internal convention that no longer exists.
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
    setdata!(mech.velocity.depthaverage_x, zero(T))
    setdata!(mech.velocity.depthaverage_y, zero(T))

    solver = PseudoTransientSolver(grid; solver_kwargs...)
    momentum isa DIVAMomentumBalance && diva_update!(mech, solver, rt, mask)

    elapsed = @elapsed result = pseudo_transient!(mech, cst_T, solver, rt, momentum, mask)

    speed = Float32.(sqrt.(
        (@views (interior(mech.velocity.depthaverage_x)[1:(end - 1), :, 1] .+
                 interior(mech.velocity.depthaverage_x)[2:end, :, 1]) ./ 2) .^ 2 .+
        (@views (interior(mech.velocity.depthaverage_y)[:, 1:(end - 1), 1] .+
                 interior(mech.velocity.depthaverage_y)[:, 2:end, 1]) ./ 2) .^ 2))

    mech = nothing
    GC.gc()
    return (; speed, elapsed, iterations = result.iterations, converged = result.converged,
           residual = result.residual)
end

const SOLVER_KWARGS = (
    abstol = 1e-3, maxiter = 2000, ncheck = 20, printout_every = 50,
    pseudo_timestep = GershgorinPseudoTimeStep(cfl = 0.99),
    convergence = ScaledResidual(),
    friction_update = ActiveFrictionUpdate(),
    tuning = AutotunedDynamicRelaxation(),
)

set_theme!(theme_latexfonts())

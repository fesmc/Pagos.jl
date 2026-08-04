#=
# Ice dynamics

Computing the field of ice velocity is the core of any ice-sheet model and also represents on of its computational bottlenecks. Therefore, many approximations have been developped, most of which are implemented in Pagos.jl via [`AbstractMomentumBalance`] and its subtypes. The choice of the dynamics solver depends on the specific research question, the available computational resources, and the desired level of accuracy. For instance, a common choice for continental ice-sheet modeling is to use a simplified dynamics solver such as the hybrid formulation of the shallow ice/shelf approximation (SIA/SSA), which is computationally efficient and can capture the large-scale flow patterns of ice sheets. However, [goldberg, robinson](@citet) have shown that the depth-integrated viscosity approximation is a more accurate and computationally efficient alternative. Therefore, it comes as the default choice in Pagos.jl, and is the recommended option for large-scale ice-sheet modeling:
=#



#=
In Pagos.jl, an important distinction is made between [`AbstractMomentumBalance`](@ref) and [`AbstractMomentumSolver`](@ref). The former represents the physical representation of the dynamics, while the latter represents the numerical implementation of the solver. This separation allows for greater flexibility and modularity in the code, as different solvers can be implemented for the same physical representation of the dynamics. For instance, one could implement a direct solver for the depth-integrated viscosity approximation, as well as an iterative solver based on a Krylov subspace method, and both would be compatible with the same physical representation of the dynamics.
=#

using Pkg
Pkg.activate(".")
using CairoMakie
using NCDatasets
using Pagos
using Statistics

absvel(ux, uy) = sqrt(ux ^ 2 + uy ^ 2)
logvel(ux, uy) = log10(absvel(ux, uy) + 1e-10)

function get_xy(X, Y; km2m = true)
    x1, x2 = extrema(X)
    y1, y2 = extrema(Y)
    lx = x2 - x1
    ly = y2 - y1
    dx = mean(diff(X, dims=1))
    dy = mean(diff(Y, dims=2))
    if km2m
        lx *= 1e3
        ly *= 1e3
        dx *= 1e3
        dy *= 1e3
    end
    return lx, ly, dx, dy
end

T = Float64
fn = "../ice-data/antarctica/topography/ANT-32KM_TOPO-RTOPO-2.0.1.nc"
ds = NCDataset(fn, "r")
mask = T.(ds["mask"][:, :])
close(ds)

fn_restart = "../ice-data/antarctica/spinups/16km/yelmo_restart.nc"
ds = NCDataset(fn_restart, "r")
H = T.(ds["H_ice"][1:2:end, 1:2:end])
zb = T.(ds["z_bed"][1:2:end, 1:2:end])
X = T.(ds["x2D"][1:2:end, 1:2:end]) .* 1e3
Y = T.(ds["y2D"][1:2:end, 1:2:end]) .* 1e3
uxy_s = T.(ds["uxy_s"][1:2:end, 1:2:end])
uxy_bar = T.(ds["uxy_bar"][1:2:end, 1:2:end])
ux_bar = T.(ds["ux_bar"][1:2:end, 1:2:end])
uy_bar = T.(ds["uy_bar"][1:2:end, 1:2:end])
cb_ref = T.(ds["cb_ref"][1:2:end, 1:2:end])
z_s = T.(ds["z_srf"][1:2:end, 1:2:end])
mu = T.(ds["visc_bar"][1:2:end, 1:2:end])
tau_b1 = T.(ds["taub_acx"][1:2:end, 1:2:end])
tau_b2 = T.(ds["taub_acy"][1:2:end, 1:2:end])
close(ds)
heatmap(H)
heatmap(z_s)
heatmap(tau_b1)
heatmap(tau_b2)
extrema(ux_bar)
extrema(uy_bar)
u_bar = cat(ux_bar, uy_bar, dims=3)

nx, ny = size(H)
dudt = zeros(T, nx, ny, 2)
u = zeros(T, nx, ny, 2)
rho_ice = 910.0
g = 9.81

t_0 = 0.0
dt = 1e-1                    # yr
t_end = 1e3               # yr
t_vec = t_0:dt:t_end
dt_save = 10.0               # yr
t_save = t_0:dt_save:t_end
u_save = zeros(T, nx, ny, 2, length(t_save))
dx, dy = 32e3, 32e3
mask = H .> 0

# inertial_velocity!(dudt, u, rho_ice, mu, H, z_s, g, tau_b1, tau_b2, dx, dy)
# @show extrema(dudt)
# u .+= dudt .* dt
# heatmap(absvel.(u[:, :, 1], u[:, :, 2]))
# @show extrema(dudt)

# inertial_velocity!(dudt, u, rho_ice, mu, H, z_s, g, tau_b1, tau_b2, dx, dy)
# u .+= dudt .* 10
# heatmap(absvel.(u[:, :, 1], u[:, :, 2]))
# @show extrema(u)
inertial_velocity!(dudt, u_bar, rho_ice, mu, H, z_s, g, tau_b1, tau_b2, dx, dy, mask)

for t in t_vec
    inertial_velocity!(dudt, u_bar, rho_ice, mu, H, z_s, g, tau_b1, tau_b2, dx, dy, mask)
    u .+= dudt .* dt
    @show extrema(u)
    if isapprox(t % dt_save, 0)
        idx = Int(t / dt_save) + 1
        u_save[:, :, :, idx] .= u
    end
end

cmap_vi = cgrad([:blue, :white, :red])
heatmap(ux_bar, colorrange = (-1e2, 1e2), colormap = cmap_vi)
heatmap(uy_bar, colorrange = (-1e2, 1e2), colormap = cmap_vi)
heatmap(absvel.(ux_bar, uy_bar), colorrange = (0, 1e3))
heatmap(u_save[:, :, 1, 20], colorrange = (-1e2, 1e2), colormap = cmap_vi)
heatmap(u_save[:, :, 2, 20], colorrange = (-1e2, 1e2), colormap = cmap_vi)
heatmap(absvel.(u[:, :, 1], u[:, :, 2]), colorrange = (0, 1e4))

II = CartesianIndices(mask)[mask]


lx, ly, dx, dy = get_xy(X, Y)
domain = Domain(T, lx, ly, dx, dy)
state = State(domain)
params = Params{T}()
options = Options{T}(
    maxiter = 10_000,
    printout_every = 1_000,
    dtau_scaling = 1e-2,
    debug = true,
)
icesheet = IceSheet(state, domain, params, options)
(;domain, state, params, options) = icesheet
mask[ 700 .< X .< 1700 .&& -1000 .< Y .< 0 ] .= 2
heatmap(mask)

state.H .= H .* (mask .== 2)
state.z_b .= zb
state.mu .= 1e5     # Pa yr
state.beta .= 1e3   # Pa yr m⁻¹
dt = 1.0            # yr

pseudo_transient!(icesheet)
logvel1 = copy(logvel.(state.ux, state.uy))
absvel1 = copy(absvel.(state.ux, state.uy))
# advect!(icesheet)
# icesheet.state.H[H .<= 0] .= 0
# logvel1_old = logvel(state.ux_old, state.uy_old)
# t2 = @elapsed pseudo_transient!(icesheet)
# logvel2 = copy(logvel(state.ux, state.uy))

nrows, ncols = 2, 3
fig = Figure(size = (1200, 900))
axs = reshape([Axis(fig[i, j], aspect = DataAspect()) for i in 1:nrows, j in 1:ncols], (nrows, ncols))
[hidedecorations!(ax) for ax in axs]

h_opts = (colormap = :ice, colorrange = (0, 4000))
dh_opts = (colormap = :balance, colorrange = (-50, 50))
# hm1 = heatmap!(axs[1, 1], H; h_opts...)
hm2 = heatmap!(axs[1, 2], state.H; h_opts...)
# hm3 = heatmap!(axs[1, 3], H - state.H; dh_opts...)
# Colorbar(fig[0, 1], hm1, vertical = false, width = Relative(0.8))
Colorbar(fig[0, 2], hm2, vertical = false, width = Relative(0.8))
# Colorbar(fig[0, 3], hm3, vertical = false, width = Relative(0.8))

log_u_opts = (colormap = :viridis, colorrange = (0, 3))
u_opts = (colormap = :viridis, colorrange = (0, 1e3))
du_opts = (colormap = :balance, colorrange = (-1, 1))
hm4 = heatmap!(axs[2, 1], logvel1; log_u_opts...)
hm5 = heatmap!(axs[2, 2], absvel1; u_opts...)
# hm6 = heatmap!(axs[2, 3], logvel2 - logvel1; du_opts...)
Colorbar(fig[3, 1], hm4, vertical = false, width = Relative(0.8))
Colorbar(fig[3, 2], hm5, vertical = false, width = Relative(0.8))
# Colorbar(fig[3, 3], hm6, vertical = false, width = Relative(0.8))

fig
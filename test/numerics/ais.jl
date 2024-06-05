using NCDatasets
using Pagos
using CairoMakie

logvel(ux, uy) = log10.(sqrt.(ux .^ 2 .+ uy .^ 2) .+ 1e-10)

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
ds = NCDataset("../ice_data/ANT-32KM_TOPO-RTOPO-2.0.1.nc", "r")
H = T.(ds["H_ice"][:, :])
zb = T.(ds["z_bed"][:, :])
X = T.(ds["x2D"][:, :])
Y = T.(ds["y2D"][:, :])
mask = T.(ds["mask"][:, :])
close(ds)

lx, ly, dx, dy = get_xy(X, Y)
domain = Domain(T, lx, ly, dx, dy)
state = State(domain)
params = Params{T}()
options = Options{T}(maxiter = 10_000, printout_every = 1_000, dtau_scaling = 1e-2, debug = true)
icesheet = IceSheet(state, domain, params, options)
(;domain, state, params, options) = icesheet
state.H .= H        # .* (mask .== 2)
state.z_b .= zb
state.mu .= 1e5     # / SEC_PER_YEAR
state.beta .= 1e3
dt = 1.0

pseudo_transient!(icesheet)
logvel1 = copy(logvel(state.ux, state.uy))
advect!(icesheet)
icesheet.state.H[H .<= 0] .= 0
# options.dtau_scaling = 1e-8
logvel1_old = logvel(state.ux_old, state.uy_old)
pseudo_transient!(icesheet)
logvel2 = copy(logvel(state.ux, state.uy))

nrows, ncols = 2, 3
fig = Figure(size = (1200, 900))
axs = reshape([Axis(fig[i, j], aspect = DataAspect()) for i in 1:nrows, j in 1:ncols], (nrows, ncols))
[hidedecorations!(ax) for ax in axs]

h_opts = (colormap = :ice, colorrange = (0, 4000))
dh_opts = (colormap = :balance, colorrange = (-50, 50))
hm1 = heatmap!(axs[1, 1], H; h_opts...)
hm2 = heatmap!(axs[1, 2], state.H; h_opts...)
hm3 = heatmap!(axs[1, 3], H - state.H; dh_opts...)
Colorbar(fig[0, 1], hm1, vertical = false, width = Relative(0.8))
Colorbar(fig[0, 2], hm2, vertical = false, width = Relative(0.8))
Colorbar(fig[0, 3], hm3, vertical = false, width = Relative(0.8))

u_opts = (colormap = :viridis, colorrange = (0, 3))
du_opts = (colormap = :balance, colorrange = (-1, 1))
hm4 = heatmap!(axs[2, 1], logvel1; u_opts...)
hm5 = heatmap!(axs[2, 2], logvel2; u_opts...)
hm6 = heatmap!(axs[2, 3], logvel2 - logvel1; du_opts...)
Colorbar(fig[3, 1], hm4, vertical = false, width = Relative(0.8))
Colorbar(fig[3, 2], hm5, vertical = false, width = Relative(0.8))
Colorbar(fig[3, 3], hm6, vertical = false, width = Relative(0.8))

fig
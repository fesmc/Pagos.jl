using Pagos
using NCDatasets
using CairoMakie

# Load topography
T = Float64
ds = NCDataset("../ice_data/yelmo2D.nc", "r")
H = T.(ds["H_ice"][:, :, :])
z_srf = T.(ds["z_srf"][:, :, :])
uxy_b = T.(ds["uxy_b"][:, :, :])
uxy_bar = T.(ds["uxy_bar"][:, :, :])
uxy_s = T.(ds["uxy_s"][:, :, :])
beta = T.(ds["beta"][:, :, :])
visc_eff_int = T.(ds["visc_eff_int"][:, :, :])
taud = T.(ds["taud"][:, :, :])
taub = T.(ds["taub"][:, :, :])
z_bed = T.(ds["z_bed"][:, :, :])
mask_ocn = T.(ds["mask_ocn"][:, :, :])
ux_s = T.(ds["ux_s"][:, :, :])
uy_s = T.(ds["uy_s"][:, :, :])
X = T.(ds["x2D"][:, :])
Y = T.(ds["y2D"][:, :])
zeta = T.(ds["zeta"][:])
t2D = T.(ds["time"][:])
close(ds)

# Missing 3D viscosity to be able to compute Fm.
ds = NCDataset("../ice_data/yelmo3D.nc", "r")
ux = T.(ds["ux"][:, :, :, :])
uy = T.(ds["uy"][:, :, :, :])
uz = T.(ds["uz"][:, :, :, :])
t3D = T.(ds["time"][:])
close(ds)

nx, ny, nz, nt3D = size(ux)

icemask = H[:, :, end] .> 0
indices2D = CartesianIndices((nx, ny))
icemask_indices = indices2D[icemask]

mu = cat([visc_eff_int[:, :, end] for l in 1:nz]..., dims=3)
b = beta[:, :, end]
sigma = zeta

# Define viscosity integrals
F1 = zeros(nx, ny)
F2 = zeros(nx, ny)
aggregated_viscosity_integral!(F1, mu, H[:, :, end], 1, sigma, icemask_indices)
aggregated_viscosity_integral!(F2, mu, H[:, :, end], 2, sigma, icemask_indices)

uxb_test = zeros(nx, ny)
uyb_test = zeros(nx, ny)
basal_velocity_from_surface_velocity!(uxb_test, uxy_s[:, :, end], b, F1, icemask_indices)
basal_velocity_from_surface_velocity!(uyb_test, uxy_s[:, :, end], b, F1, icemask_indices)

ub_test = sqrt.(uxb_test .^ 2 .+ uyb_test .^ 2)
@show extrema(ub_test) extrema(uxy_b[:, :, end])
extrema(ub_test .- uxy_b[:, :, end])
heatmap(ub_test .- uxy_b[:, :, end])
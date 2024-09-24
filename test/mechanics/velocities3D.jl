using NCDatasets
using CairoMakie
using Pagos

T = Float64
ds = NCDataset("../ice_data/ANT-32KM_TOPO-RTOPO-2.0.1.nc", "r")
H = T.(ds["H_ice"][:, :])
zb = T.(ds["z_bed"][:, :])
X = T.(ds["x2D"][:, :])
Y = T.(ds["y2D"][:, :])
mask = T.(ds["mask"][:, :])
close(ds)

ds = NCDataset("../ice_data/ANT-32KM_VEL-R11.nc", "r")
ux_srf = T.(ds["ux_srf"][:, :])
uy_srf = T.(ds["uy_srf"][:, :])
uxy_srf = T.(ds["uxy_srf"][:, :])
close(ds)

@testset "velocities3D" begin

    # Load topography
    T = Float64
    ds = NCDataset("../ice_data/ANT-32KM_TOPO-RTOPO-2.0.1.nc", "r")
    H = T.(ds["H_ice"][:, :])
    zb = T.(ds["z_bed"][:, :])
    X = T.(ds["x2D"][:, :])
    Y = T.(ds["y2D"][:, :])
    mask = T.(ds["mask"][:, :])
    close(ds)
    nx, ny = size(H)
    nz = 20
    # grounded: mask .== 2
    # floating: mask .== 3

    # Define ice mask and corresponding indices
    icemask = H .> 0
    CairoMakie.heatmap(icemask)
    indices2D = CartesianIndices((nx, ny))
    icemask_indices = indices2D[icemask]

    # Load surface velocities and crop to ice mask
    ds = NCDataset("../ice_data/ANT-32KM_VEL-R11.nc", "r")
    ux_srf = T.(ds["ux_srf"][:, :])
    uy_srf = T.(ds["uy_srf"][:, :])
    close(ds)
    ux_srf[.!icemask] .= T(0)
    uy_srf[.!icemask] .= T(0)

    # Define ice properties
    mu = fill(1e5, nx, ny, nz)              # / SEC_PER_YEAR
    beta = fill(1e3, nx, ny)
    beta[mask .== 3] .= T(0)
    sigma = exponential_vertical_layers(nz) # z3D = sigma .* H

    # Define viscosity integrals
    F1 = zeros(nx, ny)
    F2 = zeros(nx, ny)
    aggregated_viscosity_integral!(F1, mu, H, 1, sigma, icemask_indices)
    aggregated_viscosity_integral!(F2, mu, H, 2, sigma, icemask_indices)

    # Compute basal velocities
    vb_x = zeros(nx, ny)
    vb_y = zeros(nx, ny)
    vx = zeros(nx, ny)
    vy = zeros(nx, ny)
    basal_velocity_from_surface_velocity!(vb_x, ux_srf, beta, F1, icemask_indices)
    basal_velocity_from_surface_velocity!(vb_y, uy_srf, beta, F1, icemask_indices)
    vb = sqrt.(vb_x .^ 2 .+ vb_y .^ 2)
    # heatmap(vb)

    # Recover surface velocities
    vs_x = zeros(nx, ny)
    vs_y = zeros(nx, ny)
    surface_velocity!(vs_x, vb_x, beta, F1, icemask_indices)
    surface_velocity!(vs_y, vb_y, beta, F1, icemask_indices)
    @test (vs_x .≈ ux_srf) == fill(true, nx, ny)
    # heatmap(vs_x .- ux_srf)

    # Compute depth-averaged velocities
    depthavg_velocity!(vx, vb_x, beta, F2, icemask_indices)
    depthavg_velocity!(vy, vb_y, beta, F2, icemask_indices)
    v = sqrt.(vx .^ 2 .+ vy .^ 2)
    # heatmap(v)

    vs = sqrt.(ux_srf .^ 2 .+ uy_srf .^ 2)
    @test (vb .<= v) == fill(true, nx, ny)
    @test (v .<= vs) == fill(true, nx, ny)
    @show maximum(vb), maximum(v), maximum(vs)

    vx3D = zeros(nx, ny, nz)
    vy3D = zeros(nx, ny, nz)
    vz3D = zeros(nx, ny, nz)
    velocities3D!(F1, vx3D, vy3D, vz3D, vb_x, vb_y, mu, beta, H, sigma, icemask_indices)
    @show extrema(vx3D), extrema(vy3D), extrema(vz3D)
    @test (vx3D[:, :, end] .≈ ux_srf) == fill(true, nx, ny)
    @test (vy3D[:, :, end] .≈ uy_srf) == fill(true, nx, ny)

    # parabolic velocity profiles
    # lines(sigma .* H[100, 100], v3D[95, 95, :])

    # plug flow
    # lines(sigma .* H[80, 150], v3D[80, 150, :])

    # @btime velocities3D!(F1, vx3D, vy3D, vz3D, mu, beta, H, sigma)
    # 4.930 ms (160 allocations: 22.27 MiB)
    # @btime velocities3D!(F1, vx3D, vy3D, vz3D, vb_x, vb_y, mu, beta, H, sigma, icemask_indices)
    # 3.305 ms (0 allocations: 0 bytes)
end

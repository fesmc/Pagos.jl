"""
$(TYPEDSIGNATURES)

State variables for the topography component (ice geometry, mass balance, surface/bed elevations).
All boolean masks share type parameter `B`; all float fields share type parameter `M`.
"""
struct TopographyState{B, M}
    # Boolean masks
    is_ice::B
    is_ice_allowed::B
    is_grounded::B
    is_floating::B
    is_margin::B

    # Distance and grounding-zone geometry
    distance_to_margin::M
    distance_to_grline::M
    fraction_grounded::M

    # Ice and sediment thickness
    thickness_ice::M
    thickness_ice_ref::M
    thickness_ice_dt::M
    thickness_ice_effective::M
    thickness_ice_grounded::M
    thickness_sediment::M

    # Mass balance components
    massbalance_base::M
    massbalance_base_floating::M
    massbalance_base_grounded::M
    massbalance_calving_floating::M
    massbalance_calving_grounded::M
    massbalance_discharge::M
    massbalance_front::M
    massbalance_net::M
    massbalance_surface::M
    massbalance_surface_ref::M

    # Elevations
    elevation_base::M
    elevation_bed::M
    elevation_bed_ref::M
    elevation_bed_stddev::M
    elevation_seasurface::M
    elevation_surface::M
    elevation_surface_dt::M
    elevation_surface_dx::M
    elevation_surface_dy::M
end

function TopographyState(grid::RegularGrid)
    backend = KernelAbstractions.get_backend(grid.x)
    T       = eltype(grid.x)
    (; nx, ny) = grid
    b = KernelAbstractions.zeros(backend, Bool, nx, ny)
    m = KernelAbstractions.zeros(backend, T, nx, ny)
    return TopographyState(
        [copy(b) for _ in 1:5]...,   # 5 boolean fields
        [copy(m) for _ in 1:28]...,  # 28 float fields
    )
end

###############################################################

"""
$(TYPEDSIGNATURES)

State variables for the dynamics component.
`M2` is always a 2D matrix; `M23` is a 2D matrix for depth-averaged solvers (SIA, SSA)
or a 3D array for full-column solvers (Blatter–Pattyn, Stokes).
`M23` is 2D when `grid.nz == 1`, 3D otherwise.
"""
struct DynamicsState{M2, M23}
    # Depth-averaged / 2D fields
    beta::M2
    beta_eff::M2
    c_bed::M2
    v_x_bar::M2
    v_y_bar::M2
    v_x_base::M2
    v_y_base::M2
    v_x_surf::M2
    v_y_surf::M2
    v_norm_base::M2
    v_norm_surface::M2
    tau_driving::M2
    tau_base::M2
    tau_base_vertical::M2
    flux::M2
    flux_grline::M2

    # 2D or 3D fields depending on solver
    strain_effective::M23
    strain_rate_dxx::M23
    strain_rate_dyy::M23
    strain_rate_dzz::M23
    strain_rate_dxy::M23
    strain_rate_dxz::M23
    strain_rate_dyz::M23
    strain_rate_effective::M23
    stress_xx::M23
    stress_yy::M23
    stress_zz::M23
    stress_xy::M23
    stress_xz::M23
    stress_yz::M23
    stress_effective::M23
    stress_eigenvalue_1::M23
    stress_eigenvalue_2::M23
    tau_eff::M23
    tau_lateral::M23
    v_x::M23
    v_y::M23
    v_z::M23
    v_x_dx::M23
    v_x_dy::M23
    v_x_dz::M23
    v_y_dx::M23
    v_y_dy::M23
    v_y_dz::M23
    v_z_dx::M23
    v_z_dy::M23
    v_z_dz::M23
    v_norm::M23
end

function DynamicsState(grid::RegularGrid)
    backend = KernelAbstractions.get_backend(grid.x)
    T       = eltype(grid.x)
    (; nx, ny, nz) = grid
    m2  = KernelAbstractions.zeros(backend, T, nx, ny)
    m23 = nz > 1 ?
        KernelAbstractions.zeros(backend, T, nx, ny, nz) :
        KernelAbstractions.zeros(backend, T, nx, ny)
    return DynamicsState(
        [copy(m2)  for _ in 1:16]...,  # 16 depth-averaged fields
        [copy(m23) for _ in 1:32]...,  # 32 column or depth-averaged fields
    )
end

###############################################################

"""
$(TYPEDSIGNATURES)

State variables for the thermodynamics component.
`M` is a 2D matrix for column-averaged models or a 3D array for full thermodynamics.
`M` is 2D when `grid.nz == 1`, 3D otherwise.
"""
struct ThermodynamicsState{M}
    temperature_ice::M
    temperature_ice_surface::M
    temperature_ice_homologous::M
    temperature_rock::M
    temperature_pressure_melting_point::M
    enthalpy_ice::M
    enthalpy_rock::M
    ice_water_content::M
    heat_strain_internal::M
    heat_strain_internal_dt::M
    heat_base_friction::M
    heatflux_base_ice::M
    heatflux_bedrock::M
    heatflux_geothermal::M
    thickness_waterlayer::M
    thickness_waterlayer_dt::M
    thickness_coldtemperate_interface::M
    specific_heat_capacity_ice::M
    heat_conductivity_ice::M
end

function ThermodynamicsState(grid::RegularGrid)
    backend = KernelAbstractions.get_backend(grid.x)
    T       = eltype(grid.x)
    (; nx, ny, nz) = grid
    m = nz > 1 ?
        KernelAbstractions.zeros(backend, T, nx, ny, nz) :
        KernelAbstractions.zeros(backend, T, nx, ny)
    return ThermodynamicsState([copy(m) for _ in 1:19]...)
end

###############################################################

"""
$(TYPEDSIGNATURES)

State variables for the material component.
`M2` is a depth-averaged 2D matrix; `M3` is a full 3D array.
"""
struct MaterialState{M2, M3}
    eta_depth_averaged::M2
    eta_depth_integrated::M2
    eta_ice::M3
end

function MaterialState(grid::RegularGrid)
    backend = KernelAbstractions.get_backend(grid.x)
    T       = eltype(grid.x)
    (; nx, ny, nz) = grid
    m2 = KernelAbstractions.zeros(backend, T, nx, ny)
    m3 = KernelAbstractions.zeros(backend, T, nx, ny, nz)
    return MaterialState(copy(m2), copy(m2), m3)
end

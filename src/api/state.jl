struct TopographicMasks{B}
    is_ice::B
    is_ice_allowed::B
    is_grounded::B
    is_floating::B
    is_margin::B
end
Adapt.@adapt_structure TopographicMasks

struct DistanceState{M}
    distance_to_margin::M
    distance_to_grline::M
end
Adapt.@adapt_structure DistanceState

struct FractionState{M}
    fraction_grounded::M
end
Adapt.@adapt_structure FractionState

struct MassBalanceState{M}
    base::M
    base_floating::M
    base_grounded::M
    calving_floating::M
    calving_grounded::M
    discharge::M
    front::M
    net::M
    surface::M
    surface_ref::M
end
Adapt.@adapt_structure MassBalanceState

struct ElevationState{M}
    base::M
    bed::M
    bed_ref::M
    bed_stddev::M
    seasurface::M
    surface::M
    surface_dt::M
    surface_dx::M
    surface_dy::M
end
Adapt.@adapt_structure ElevationState

struct ThicknessState{M}
    ice::M
    ice_ref::M
    ice_dt::M
    ice_effective::M
    ice_grounded::M
    sediment::M
end
Adapt.@adapt_structure ThicknessState

"""
$(TYPEDSIGNATURES)

State variables for the topography component (ice geometry, mass balance, surface/bed elevations).
All boolean masks share type parameter `B`; all float fields share type parameter `M`.
"""
struct TopographicState{B, M}
    mask::TopographicMasks{B}
    distance::DistanceState{M}
    fraction::FractionState{M}
    thickness::ThicknessState{M}
    massbalance::MassBalanceState{M}
    elevation::ElevationState{M}
end
Adapt.@adapt_structure TopographicState


function TopographicState(grid::RegularGrid)
    backend = KernelAbstractions.get_backend(grid.x)
    T       = eltype(grid.x)
    (; nx, ny) = grid
    b = KernelAbstractions.zeros(backend, Bool, nx, ny)
    m = KernelAbstractions.zeros(backend, T, nx, ny)
    return TopographicState(
        TopographicMasks([copy(b) for _ in 1:5]...),
        DistanceState([copy(m) for _ in 1:2]...),
        FractionState([copy(m) for _ in 1:1]...),
        ThicknessState([copy(m) for _ in 1:6]...),
        MassBalanceState([copy(m) for _ in 1:10]...),
        ElevationState([copy(m) for _ in 1:9]...),
    )
end

###############################################################

struct MechanicTopographyState{M}
    surface::M
    thickness::M
end
Adapt.@adapt_structure MechanicTopographyState

struct MechanicMaterialState{M2, M3}
    viscosity_depthaveraged::M2
    viscosity::M3
end
Adapt.@adapt_structure MechanicMaterialState

struct StressState{M2, M3}
    driving::M2
    base::M2
    base_vertical::M2

    xx::M3
    xy::M3
    xz::M3
    yx::M3
    yy::M3
    yz::M3
    zx::M3
    zy::M3
    zz::M3
    effective::M3
    lateral::M3
    eigenvalue_1::M3
    eigenvalue_2::M3
end
Adapt.@adapt_structure StressState

struct StrainRateState{M3}
    xx::M3
    xy::M3
    xz::M3
    yx::M3
    yy::M3
    yz::M3
    zx::M3
    zy::M3
    zz::M3
    effective::M3
end
Adapt.@adapt_structure StrainRateState

struct VelocityState{M2, M3}
    x_bar::M2
    y_bar::M2
    x_bar_dx::M2
    x_bar_dy::M2
    x_bar_dz::M2
    y_bar_dx::M2
    y_bar_dy::M2
    y_bar_dz::M2

    x_base::M2
    y_base::M2
    x_surf::M2
    y_surf::M2
    norm_base::M2
    norm_surface::M2

    x::M3
    y::M3
    z::M3
    x_dx::M3
    x_dy::M3
    x_dz::M3
    y_dx::M3
    y_dy::M3
    y_dz::M3
    z_dx::M3
    z_dy::M3
    z_dz::M3
    norm::M3
end
Adapt.@adapt_structure VelocityState

"""
$(TYPEDSIGNATURES)

State variables for the dynamics component.
`M2` is always a 2D matrix `(nx, ny)` for the depth-averaged / vertically-integrated
fields. `M3` is always a 3D array `(nx, ny, nz)` for the column fields (velocity, its
gradients, the strain-rate and stress tensors). Depth-averaged solvers (SIA, SSA) simply
use `nz == 1`, so the column dimension is always present and tensor/stress computations
are dynamics-independent (no 2D/3D special-casing).
"""
struct MechanicState{M2, M3}
    beta::M2
    beta_eff::M2
    c_bed::M2
    flux::M2
    flux_grline::M2

    # Column (3D) fields; depth-averaged solvers use nz == 1
    topography::MechanicTopographyState{M2}
    material::MechanicMaterialState{M2, M3}
    strainrate::StrainRateState{M3}
    stress::StressState{M2, M3}
    velocity::VelocityState{M2, M3}
end
Adapt.@adapt_structure MechanicState

function MechanicState(grid::RegularGrid)
    backend = KernelAbstractions.get_backend(grid.x)
    T       = eltype(grid.x)
    (; nx, ny, nz) = grid
    m2() = KernelAbstractions.zeros(backend, T, nx, ny)
    m3() = KernelAbstractions.zeros(backend, T, nx, ny, nz)
    return MechanicState(
        m2(), m2(), m2(), m2(), m2(),                       # beta, beta_eff, c_bed, flux, flux_grline
        MechanicTopographyState(m2(), m2()),               # surface, thickness
        MechanicMaterialState(m2(), m3()),                 # viscosity_depthaveraged, viscosity
        StrainRateState(ntuple(_ -> m3(), 10)...),         # 10 column tensor fields
        StressState(m2(), m2(), m2(), ntuple(_ -> m3(), 13)...),  # 3 depth-averaged + 13 column
        VelocityState(ntuple(_ -> m2(), 14)..., ntuple(_ -> m3(), 13)...),  # 14 depth-averaged + 13 column
    )
end

###############################################################

struct TemperatureState{M}
    ice::M
    ice_surface::M
    ice_homologous::M
    rock::M
    pressure_melting_point::M
end
Adapt.@adapt_structure TemperatureState

struct EnthalpyState{M}
    ice::M
    rock::M
end
Adapt.@adapt_structure EnthalpyState

"""
$(TYPEDSIGNATURES)

State variables for the thermodynamics component.
`M` is a 2D matrix for column-averaged models or a 3D array for full thermodynamics.
`M` is 2D when `grid.nz == 1`, 3D otherwise.
"""
struct ThermodynamicState{M}
    temperature::TemperatureState{M}
    enthalpy::EnthalpyState{M}
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
Adapt.@adapt_structure ThermodynamicState

function ThermodynamicState(grid::RegularGrid)
    backend = KernelAbstractions.get_backend(grid.x)
    T       = eltype(grid.x)
    (; nx, ny, nz) = grid
    m() = nz > 1 ?
        KernelAbstractions.zeros(backend, T, nx, ny, nz) :
        KernelAbstractions.zeros(backend, T, nx, ny)
    return ThermodynamicState(
        TemperatureState(ntuple(_ -> m(), 5)...),  # ice, ice_surface, ice_homologous, rock, pressure_melting_point
        EnthalpyState(m(), m()),                   # ice, rock
        ntuple(_ -> m(), 12)...,                   # 12 remaining flat fields
    )
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
Adapt.@adapt_structure MaterialState

function MaterialState(grid::RegularGrid)
    backend = KernelAbstractions.get_backend(grid.x)
    T       = eltype(grid.x)
    (; nx, ny, nz) = grid
    m2 = KernelAbstractions.zeros(backend, T, nx, ny)
    m3 = KernelAbstractions.zeros(backend, T, nx, ny, nz)
    return MaterialState(copy(m2), copy(m2), m3)
end

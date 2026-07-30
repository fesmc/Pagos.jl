###############################################################
# Field layout: the Arakawa C-grid node classes
###############################################################
#
# Chmy expresses staggering as a per-axis location, so the Yelmo node classes the field
# names already gesture at (`taud_acx`, `β_acx`, ...) become location tuples:
#
#   | Yelmo node | horizontal location  | typical fields                              |
#   |------------|----------------------|---------------------------------------------|
#   | `aa`       | `(Center, Center)`   | H, z_srf, z_bed, T, A, η, ε̇_xx, ε̇_yy       |
#   | `acx`      | `(Vertex, Center)`   | u, taud_acx, β_acx, q_x                     |
#   | `acy`      | `(Center, Vertex)`   | v, taud_acy, β_acy, q_y                     |
#   | `ab`       | `(Vertex, Vertex)`   | ε̇_xy, σ_xy, corner viscosity                |
#
# Vertically, layer midpoints (`ζ_aa`) are z-`Center` and layer interfaces (`ζ_ac`) are
# z-`Vertex`; the `_AC` suffix below marks the interface variants (vertical shear ε̇_xz,
# ε̇_yz, the vertical velocity w and everything differentiated with respect to z).
#
# The locations are not free choices: they follow from the operators. `∂x` maps
# `Vertex → Center` along x and leaves the other axes alone, so e.g. `∂u/∂y` with
# `u` at `acx` lands on `ab`, exactly where ε̇_xy lives, and `∂w/∂x` with `w` at
# `aa_ac` lands on `acx_ac`, exactly where ε̇_xz lives. The assignment below is the
# unique one consistent with the C-grid (and agrees with Chmy's own `TensorField{3}`
# component locations).

const NODE_AA      = (Center(), Center(), Center())
const NODE_ACX     = (Vertex(), Center(), Center())
const NODE_ACY     = (Center(), Vertex(), Center())
const NODE_AB      = (Vertex(), Vertex(), Center())
const NODE_AA_AC   = (Center(), Center(), Vertex())
const NODE_ACX_AC  = (Vertex(), Center(), Vertex())
const NODE_ACY_AC  = (Center(), Vertex(), Vertex())

"""
$(TYPEDSIGNATURES)

Allocate a Chmy `Field` of element type `T` at location `loc` on `grid`, for the
architecture `arch`. Thin, eltype-explicit wrapper around `Chmy.Field` used by the
`StaggeredGrid` methods of the state constructors.
"""
_field(arch, grid, loc, T, halo) = Field(arch, grid, loc, T; halo)

###############################################################

struct TopographicMasks{AA}
    is_ice::AA
    is_ice_allowed::AA
    is_grounded::AA
    is_floating::AA
    is_margin::AA
end
Adapt.@adapt_structure TopographicMasks

struct DistanceState{AA}
    distance_to_margin::AA
    distance_to_grline::AA
end
Adapt.@adapt_structure DistanceState

struct FractionState{AA}
    fraction_grounded::AA
end
Adapt.@adapt_structure FractionState

struct MassBalanceState{AA}
    base::AA
    base_floating::AA
    base_grounded::AA
    calving_floating::AA
    calving_grounded::AA
    discharge::AA
    front::AA
    net::AA
    surface::AA
    surface_ref::AA
end
Adapt.@adapt_structure MassBalanceState

struct ThicknessState{AA}
    ice::AA
    ice_ref::AA
    ice_dt::AA
    ice_effective::AA
    ice_grounded::AA
    sediment::AA
end
Adapt.@adapt_structure ThicknessState

struct ElevationState{AA, ACX, ACY}
    base::AA
    bed::AA
    bed_ref::AA
    bed_stddev::AA
    seasurface::AA
    surface::AA
    surface_dt::AA
    surface_dx::ACX     # ∂s/∂x lives on the acx face, where the driving stress does
    surface_dy::ACY     # ∂s/∂y lives on the acy face
end
Adapt.@adapt_structure ElevationState

"""
$(TYPEDSIGNATURES)

State variables for the topography component (ice geometry, mass balance, surface/bed
elevations). Every field is depth-integrated: on a [`StaggeredGrid`](@ref) they are all
built on `grid.grid2d`.

Type parameters are array (or `Chmy.Field`) types named after the node class they carry:
`B` are the boolean masks at `aa`, `M` the float fields at `aa`, `ACX`/`ACY` the surface
gradients on the `acx`/`acy` faces. With plain arrays ([`RegularGrid`](@ref)) `M`, `ACX`
and `ACY` are the same type; with `Field`s they differ, because the location is part of
the type.
"""
struct TopographicState{B, M, ACX, ACY}
    mask::TopographicMasks{B}
    distance::DistanceState{M}
    fraction::FractionState{M}
    thickness::ThicknessState{M}
    massbalance::MassBalanceState{M}
    elevation::ElevationState{M, ACX, ACY}
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

"""
$(TYPEDSIGNATURES)

Build a [`TopographicState`](@ref) of Chmy `Field`s on `grid`. All fields are
depth-integrated, so they live on `grid.grid2d` at the `aa` node — except the surface
gradients, which live on the `acx`/`acy` faces. Element type follows the grid; `halo`
is the ghost-cell width of every field.
"""
function TopographicState(grid::StaggeredGrid; halo = 1)
    (; arch) = grid
    g = grid.grid2d
    T = eltype(g)
    b()   = _field(arch, g, NODE_AA, Bool, halo)
    aa()  = _field(arch, g, NODE_AA, T, halo)
    acx() = _field(arch, g, NODE_ACX, T, halo)
    acy() = _field(arch, g, NODE_ACY, T, halo)
    return TopographicState(
        TopographicMasks(ntuple(_ -> b(), 5)...),
        DistanceState(ntuple(_ -> aa(), 2)...),
        FractionState(aa()),
        ThicknessState(ntuple(_ -> aa(), 6)...),
        MassBalanceState(ntuple(_ -> aa(), 10)...),
        ElevationState(ntuple(_ -> aa(), 7)..., acx(), acy()),
    )
end

###############################################################

struct MechanicTopographyState{AA}
    surface::AA
    thickness::AA
end
Adapt.@adapt_structure MechanicTopographyState

struct MechanicMaterialState{AA2, AA3}
    viscosity_depthaveraged::AA2
    viscosity::AA3
end
Adapt.@adapt_structure MechanicMaterialState

# `beta`, `beta_eff` and `c_bed` are evaluated at `aa`, where the effective pressure and
# the basal velocity magnitude they depend on live. The C-grid needs β on the velocity
# faces (`β_acx`, `β_acy`); whether those become stored fields or an inline `lerp` inside
# the momentum kernel is a Phase 3 decision (see `roadmaps/chmy.md`).
struct FrictionState{AA}
    beta::AA
    beta_eff::AA
    c_bed::AA
end
Adapt.@adapt_structure FrictionState

# The depth-integrated mass flux q = H ū is what the continuity equation differentiates,
# and `∂H/∂t = -divg(q)` is only conservative if the two components sit on the cell faces
# the divergence reads (`acx`/`acy`) — a single `aa` flux would have to be re-staggered
# inside the divergence, which is exactly the half-cell shift the C-grid exists to avoid.
# `grline` is a grounding-line diagnostic, not a term in the continuity equation, so it
# stays a cell-centred scalar at `aa`.
struct FluxState{ACX2, ACY2, AA2}
    x::ACX2
    y::ACY2
    grline::AA2
end
Adapt.@adapt_structure FluxState

struct StressState{ACX2, ACY2, AA2, AA3, AB3, ACXZ3, ACYZ3}
    driving_x::ACX2
    driving_y::ACY2
    base_x::ACX2
    base_y::ACY2
    base_vertical::AA2

    xx::AA3
    xy::AB3
    xz::ACXZ3
    yx::AB3
    yy::AA3
    yz::ACYZ3
    zx::ACXZ3
    zy::ACYZ3
    zz::AA3
    effective::AA3
    lateral::AA3
    eigenvalue_1::AA3
    eigenvalue_2::AA3
end
Adapt.@adapt_structure StressState

struct StrainRateState{AA3, AB3, ACXZ3, ACYZ3}
    xx::AA3
    xy::AB3
    xz::ACXZ3
    yx::AB3
    yy::AA3
    yz::ACYZ3
    zx::ACXZ3
    zy::ACYZ3
    zz::AA3
    effective::AA3
end
Adapt.@adapt_structure StrainRateState

struct VelocityState{ACX2, ACY2, AA2, AB2, ACX3, ACY3, AA3, AB3, AAZ3, ACXZ3, ACYZ3}
    x_bar::ACX2
    y_bar::ACY2
    x_bar_dx::AA2
    x_bar_dy::AB2
    x_bar_dz::ACX2
    y_bar_dx::AB2
    y_bar_dy::AA2
    y_bar_dz::ACY2

    x_base::ACX2
    y_base::ACY2
    x_surf::ACX2
    y_surf::ACY2
    norm_base::AA2
    norm_surface::AA2

    x::ACX3
    y::ACY3
    z::AAZ3
    x_dx::AA3
    x_dy::AB3
    x_dz::ACXZ3
    y_dx::AB3
    y_dy::AA3
    y_dz::ACYZ3
    z_dx::ACXZ3
    z_dy::ACYZ3
    z_dz::AA3
    norm::AA3
end
Adapt.@adapt_structure VelocityState

"""
$(TYPEDSIGNATURES)

State variables for the dynamics component.

Type parameters are the array (or `Chmy.Field`) types of the C-grid node classes, `2`
marking a depth-integrated field and `3` a column field: `AA2`/`AA3` at cell centres,
`ACX*`/`ACY*` on the x-/y-faces, `AB*` on the corners, `AAZ3`/`ACXZ3`/`ACYZ3` on the
layer interfaces (`ζ_ac`) at the corresponding horizontal node. With plain arrays
([`RegularGrid`](@ref)) all the `*2` parameters collapse to one matrix type and all the
`*3` parameters to one 3-array type; with `Field`s they are distinct types, so a
wrong-location assignment is a construction error rather than a silent half-cell shift.

Column fields carry the vertical dimension even for depth-averaged momentum balances
(SIA, SSA), which simply use `nz == 1`, so tensor and stress computations stay
dynamics-independent (no 2D/3D special-casing).
"""
struct MechanicState{ACX2, ACY2, AA2, AB2, ACX3, ACY3, AA3, AB3, AAZ3, ACXZ3, ACYZ3}
    friction::FrictionState{AA2}
    flux::FluxState{ACX2, ACY2, AA2}

    topography::MechanicTopographyState{AA2}
    material::MechanicMaterialState{AA2, AA3}
    strainrate::StrainRateState{AA3, AB3, ACXZ3, ACYZ3}
    stress::StressState{ACX2, ACY2, AA2, AA3, AB3, ACXZ3, ACYZ3}
    velocity::VelocityState{ACX2, ACY2, AA2, AB2, ACX3, ACY3, AA3, AB3, AAZ3, ACXZ3, ACYZ3}
end
Adapt.@adapt_structure MechanicState

function MechanicState(grid::RegularGrid)
    backend = KernelAbstractions.get_backend(grid.x)
    T       = eltype(grid.x)
    (; nx, ny, nz) = grid
    m2() = KernelAbstractions.zeros(backend, T, nx, ny)
    m3() = KernelAbstractions.zeros(backend, T, nx, ny, nz)
    return MechanicState(
        FrictionState(m2(), m2(), m2()),                    # beta, beta_eff, c_bed
        FluxState(m2(), m2(), m2()),                       # flux x, y, grline
        MechanicTopographyState(m2(), m2()),               # surface, thickness
        MechanicMaterialState(m2(), m3()),                 # viscosity_depthaveraged, viscosity
        StrainRateState(ntuple(_ -> m3(), 10)...),         # 10 column tensor fields
        StressState(ntuple(_ -> m2(), 5)..., ntuple(_ -> m3(), 13)...),  # 3 depth-averaged + 13 column
        VelocityState(ntuple(_ -> m2(), 14)..., ntuple(_ -> m3(), 13)...),  # 14 depth-averaged + 13 column
    )
end

"""
$(TYPEDSIGNATURES)

Build a [`MechanicState`](@ref) of Chmy `Field`s on `grid`. Depth-integrated fields are
built on `grid.grid2d`, column fields on `grid.grid`; each is placed at the node class
its physics dictates (see the table at the top of `src/api/state.jl`). Element type
follows the grid; `halo` is the ghost-cell width of every field.
"""
function MechanicState(grid::StaggeredGrid; halo = 1)
    (; arch) = grid
    g2, g3 = grid.grid2d, grid.grid
    T = eltype(g3)
    aa2()  = _field(arch, g2, NODE_AA, T, halo)
    acx2() = _field(arch, g2, NODE_ACX, T, halo)
    acy2() = _field(arch, g2, NODE_ACY, T, halo)
    ab2()  = _field(arch, g2, NODE_AB, T, halo)
    aa3()  = _field(arch, g3, NODE_AA, T, halo)
    acx3() = _field(arch, g3, NODE_ACX, T, halo)
    acy3() = _field(arch, g3, NODE_ACY, T, halo)
    ab3()  = _field(arch, g3, NODE_AB, T, halo)
    aaz3()  = _field(arch, g3, NODE_AA_AC, T, halo)
    acxz3() = _field(arch, g3, NODE_ACX_AC, T, halo)
    acyz3() = _field(arch, g3, NODE_ACY_AC, T, halo)
    return MechanicState(
        FrictionState(aa2(), aa2(), aa2()),
        FluxState(acx2(), acy2(), aa2()),                  # flux x, y, grline
        MechanicTopographyState(aa2(), aa2()),
        MechanicMaterialState(aa2(), aa3()),
        StrainRateState(
            aa3(), ab3(), acxz3(), ab3(), aa3(),         # xx, xy, xz, yx, yy
            acyz3(), acxz3(), acyz3(), aa3(), aa3(),     # yz, zx, zy, zz, effective
        ),
        StressState(
            acx2(), acy2(), acx2(), acy2(), aa2(),       # driving_x/y, base_x/y, base_vertical
            aa3(), ab3(), acxz3(), ab3(), aa3(),         # xx, xy, xz, yx, yy
            acyz3(), acxz3(), acyz3(), aa3(),            # yz, zx, zy, zz
            aa3(), aa3(), aa3(), aa3(),                  # effective, lateral, eigenvalue_1/2
        ),
        VelocityState(
            acx2(), acy2(), aa2(), ab2(), acx2(),        # x_bar, y_bar, x_bar_dx/dy/dz
            ab2(), aa2(), acy2(),                        # y_bar_dx/dy/dz
            acx2(), acy2(), acx2(), acy2(), aa2(), aa2(),  # x/y_base, x/y_surf, norm_base/surface
            acx3(), acy3(), aaz3(),                      # x, y, z
            aa3(), ab3(), acxz3(),                       # x_dx, x_dy, x_dz
            ab3(), aa3(), acyz3(),                       # y_dx, y_dy, y_dz
            acxz3(), acyz3(), aa3(), aa3(),              # z_dx, z_dy, z_dz, norm
        ),
    )
end

###############################################################

struct TemperatureState{AA3, AA2}
    ice::AA3
    ice_surface::AA2
    ice_homologous::AA3
    rock::AA3
    pressure_melting_point::AA3
end
Adapt.@adapt_structure TemperatureState

struct EnthalpyState{AA3}
    ice::AA3
    rock::AA3
end
Adapt.@adapt_structure EnthalpyState

"""
$(TYPEDSIGNATURES)

State variables for the thermodynamics component. `AA3` are the column fields at `aa`
(temperature, enthalpy, water content, the internal strain heating and the thermal
material properties), `AA2` the depth-integrated ones — surface temperature, the heat
fluxes at the ice base and in the bedrock, the basal water layer and the position of the
cold-temperate transition surface, none of which carry a vertical index.

On a [`RegularGrid`](@ref) both parameters collapse to a single array type. On a
[`StaggeredGrid`](@ref) the `AA3` fields are built on `grid.grid` and the `AA2` fields on
`grid.grid2d`.
"""
struct ThermodynamicState{AA3, AA2}
    temperature::TemperatureState{AA3, AA2}
    enthalpy::EnthalpyState{AA3}
    ice_water_content::AA3
    heat_strain_internal::AA3
    heat_strain_internal_dt::AA3
    heat_base_friction::AA2
    heatflux_base_ice::AA2
    heatflux_bedrock::AA2
    heatflux_geothermal::AA2
    thickness_waterlayer::AA2
    thickness_waterlayer_dt::AA2
    thickness_coldtemperate_interface::AA2
    specific_heat_capacity_ice::AA3
    heat_conductivity_ice::AA3
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

"""
$(TYPEDSIGNATURES)

Build a [`ThermodynamicState`](@ref) of Chmy `Field`s on `grid`, all at the `aa` node:
column fields on `grid.grid`, depth-integrated ones on `grid.grid2d` (see the
[`ThermodynamicState`](@ref) docstring for which is which). Element type follows the
grid; `halo` is the ghost-cell width of every field.
"""
function ThermodynamicState(grid::StaggeredGrid; halo = 1)
    (; arch) = grid
    g2, g3 = grid.grid2d, grid.grid
    T = eltype(g3)
    aa2() = _field(arch, g2, NODE_AA, T, halo)
    aa3() = _field(arch, g3, NODE_AA, T, halo)
    return ThermodynamicState(
        TemperatureState(aa3(), aa2(), aa3(), aa3(), aa3()),
        EnthalpyState(aa3(), aa3()),
        aa3(), aa3(), aa3(),                       # ice_water_content, heat_strain_internal(_dt)
        aa2(), aa2(), aa2(), aa2(),                # heat_base_friction, heatflux_base_ice/bedrock/geothermal
        aa2(), aa2(), aa2(),                       # thickness_waterlayer(_dt), thickness_coldtemperate_interface
        aa3(), aa3(),                              # specific_heat_capacity_ice, heat_conductivity_ice
    )
end

###############################################################

"""
$(TYPEDSIGNATURES)

State variables for the material component, all at the `aa` node: `AA2` are the
depth-averaged and depth-integrated viscosities, `AA3` the viscosity of the ice column.
"""
struct MaterialState{AA2, AA3}
    eta_depth_averaged::AA2
    eta_depth_integrated::AA2
    eta_ice::AA3
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

"""
$(TYPEDSIGNATURES)

Build a [`MaterialState`](@ref) of Chmy `Field`s on `grid`: the depth-averaged and
depth-integrated viscosities on `grid.grid2d`, the column viscosity on `grid.grid`, all
at the `aa` node. Element type follows the grid; `halo` is the ghost-cell width of every
field.
"""
function MaterialState(grid::StaggeredGrid; halo = 1)
    (; arch) = grid
    T = eltype(grid.grid)
    aa2() = _field(arch, grid.grid2d, NODE_AA, T, halo)
    return MaterialState(aa2(), aa2(), _field(arch, grid.grid, NODE_AA, T, halo))
end

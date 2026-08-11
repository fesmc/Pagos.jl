###############################################################
# Field layout: the Arakawa C-grid node classes
###############################################################
# Chmy expresses staggering as a per-axis location; Yelmo's node names map to location
# tuples:
#
#   | Yelmo node | horizontal location  | typical fields                              |
#   |------------|----------------------|---------------------------------------------|
#   | `aa`       | `(Center, Center)`   | H, z_srf, z_bed, T, A, η, ε̇_xx, ε̇_yy      |
#   | `acx`      | `(Vertex, Center)`   | u, taud_acx, β_acx, q_x                     |
#   | `acy`      | `(Center, Vertex)`   | v, taud_acy, β_acy, q_y                     |
#   | `ab`       | `(Vertex, Vertex)`   | ε̇_xy, σ_xy, corner viscosity               |
#
# Vertically, layer midpoints (`ζ_aa`) are z-`Center`, layer interfaces (`ζ_ac`) are
# z-`Vertex`; the `_AC` suffix below marks the interface variants (vertical shear ε̇_xz,
# ε̇_yz, w, and anything differentiated w.r.t. z).
#
# Locations follow from the operators, not free choice: `∂x` maps `Vertex → Center` along
# x, so `∂u/∂y` with `u` at `acx` lands on `ab` (where ε̇_xy lives), and `∂w/∂x` with `w`
# at `aa_ac` lands on `acx_ac` (where ε̇_xz lives) — the unique assignment consistent with
# the C-grid, matching Chmy's own `TensorField{3}` component locations.

const NODE_AA = (Center(), Center(), Center())
const NODE_ACX = (Vertex(), Center(), Center())
const NODE_ACY = (Center(), Vertex(), Center())
const NODE_AB = (Vertex(), Vertex(), Center())
const NODE_AA_AC = (Center(), Center(), Vertex())
const NODE_ACX_AC = (Vertex(), Center(), Vertex())
const NODE_ACY_AC = (Center(), Vertex(), Vertex())

"""
$(TYPEDSIGNATURES)

Allocate a Chmy `Field` of element type `T` at location `loc` on `grid`, for the
architecture `arch`. Thin, eltype-explicit wrapper around `Chmy.Field` used by the
`StaggeredGrid` methods of the state constructors.
"""
_field(arch, grid, loc, T, halo) = Field(arch, grid, loc, T; halo)

"""
$(TYPEDSIGNATURES)

The halo width to allocate a field on `grid.grid2d` with: whatever the caller asked for on
`x`/`y`, and **zero** on `z`.

`grid2d`'s z-axis has extent 1 and nothing sweeps a ghost ring there, but `Chmy.Field` does
not know that: a scalar `halo` allocates the usual `4·halo` extra z-planes anyway, **5.1× the
memory the interior needs** at the default `halo = 1` (`pagos-roadmap/chmy-issues.md`). Every
`grid2d`-based field builder below passes this instead of `halo` directly, so a given `halo`
always yields the same `Field{T,N,L,H,A}` type.
"""
_halo2d(halo::Integer) = (halo, halo, 0)
_halo2d(halo::NTuple{3,Integer}) = (halo[1], halo[2], 0)

###############################################################
# Diagnostic (output-only) fields
###############################################################
#
# Fields marked `@diagnostic` below are those `docs/src/variables.md` classifies as neither
# Necessary nor Used: nothing in `src/` touches them, they exist as a place to put a quantity
# somebody may want to *output*. At a 1522² Antarctic setup they cost ~1.4 GiB across
# `TopographicState` + `MechanicState`, so `diagnostics = false` (the default) allocates them
# on a one-cell grid instead. Pass `diagnostics = true` for full size.
#
# The struct field itself stays: a `Field` over a 1-cell grid has the *same* `Field{T,N,L,H,A}`
# type as one over the full grid, so the single-type-parameter states keep working and
# `variables.md` stays an accurate description of the design.

_topology_type(::StructuredGrid{N,T,C}) where {N,T,C} = C

"""
$(TYPEDSIGNATURES)

A one-cell `StaggeredGrid` on the same architecture and element type as `grid`, used to give
[diagnostic fields](@ref variables) a real `Field` of the correct type at negligible cost.

Mirrors `StaggeredGrid`'s own axis/topology construction so the resulting fields are the same
Julia type as their full-size counterparts — only smaller.
"""
function _degenerate_grid(grid::StaggeredGrid)
    (; arch) = grid
    T = eltype(grid.grid)
    ax = UniformAxis(T(0), T(1), 1)
    return StructuredGrid{_topology_type(grid.grid)}(arch, ax, ax, ax)
end

# `_field` for a field that is diagnostic-only: full size when `diagnostics`, else one cell.
_maybe_field(arch, grid, degenerate, loc, T, halo, diagnostics) =
    _field(arch, diagnostics ? grid : degenerate, loc, T, halo)

###############################################################

"""
$(TYPEDSIGNATURES)

## Fields
- `is_ice`: ice-covered cells
- `is_ice_neighbour`: ice-free but touching ice: the ring the margin advances into
- `is_ice_allowed`: ice is allowed to exist here
- `is_grounded`: ice is grounded here
- `is_floating`: ice is floating here
- `is_margin`: cell has at least one ice-free neighbour
- `is_momentum_solved`: the momentum balance is well-posed here — ice connected to
  grounded ice through ice, i.e. everything except detached icebergs. Written by
  [`momentum_mask!`](@ref); see its docstring for why this is not the same as `is_ice`.
"""
struct TopographyMasks{AA}
    is_ice::AA
    is_ice_neighbour::AA
    is_ice_allowed::AA
    is_grounded::AA
    is_floating::AA
    is_margin::AA
    is_momentum_solved::AA
end
Adapt.@adapt_structure TopographyMasks

"""
$(TYPEDSIGNATURES)

Allocate a standalone [`TopographyMasks`](@ref) of Chmy `Field`s on `grid`, at the `aa`
node. For callers that only need the ice/grounding masks — e.g. to run [`icemasks!`](@ref)
and [`momentum_mask!`](@ref) — without paying for the rest of a
[`TopographicState`](@ref).
"""
function TopographyMasks(grid::StaggeredGrid; halo = 1)
    (; arch) = grid
    g = grid.grid2d
    return TopographyMasks(ntuple(_ -> _field(arch, g, NODE_AA, Bool, _halo2d(halo)), 7)...)
end

"""
$(TYPEDSIGNATURES)

## Fields
- `distance_to_margin`: distance to the nearest ice-free cell
- `distance_to_grline`: distance to the nearest grounding line cell
"""
struct DistanceState{AA}
    distance_to_margin::AA
    distance_to_grline::AA
end
Adapt.@adapt_structure DistanceState

"""
$(TYPEDSIGNATURES)

## Fields
- `fraction_grounded`: fraction of the cell that is grounded ice
"""
struct FractionState{AA}
    fraction_grounded::AA
end
Adapt.@adapt_structure FractionState

"""
$(TYPEDSIGNATURES)

## Fields
- `base`: basal mass balance (m/yr)
- `base_floating`: basal mass balance from subshelf melting (m/yr)
- `base_grounded`: basal mass balance from strain heating, geothermal flux, etc. (m/yr)
- `calving_floating`: calving mass loss from floating ice (m/yr)
- `calving_grounded`: calving mass loss from grounded ice (m/yr)
- `discharge`: mass loss from ice flowing out of the grounded domain (m/yr)
- `front`: mass loss from ice flowing out of the floating domain (m/yr)
- `net`: net mass balance (m/yr)
- `surface`: surface mass balance (m/yr)
- `surface_ref`: reference surface mass balance (m/yr)
"""
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

"""
$(TYPEDSIGNATURES)

## Fields
- `ice`: ice thickness (m)
- `ice_ref`: reference ice thickness (m)
- `ice_dt`: ice thickness change rate (m/yr)
- `ice_effective`: effective ice thickness (m)
- `ice_grounded`: grounded ice thickness (m)
- `sediment`: sediment thickness (m)
"""
struct ThicknessState{AA}
    ice::AA
    ice_ref::AA
    ice_dt::AA
    ice_effective::AA
    ice_grounded::AA
    sediment::AA
end
Adapt.@adapt_structure ThicknessState

"""
$(TYPEDSIGNATURES)

## Fields

- `base`: ice base elevation (m); differs from bed elevation when sediments are present or when ice is floating
- `bed`: bedrock elevation (m)
- `bed_ref`: reference bedrock elevation (m)
- `bed_stddev`: bedrock elevation standard deviation (m)
- `seasurface`: sea surface elevation (m)
- `surface`: ice surface elevation (m)
- `surface_dt`: ice surface elevation change rate (m/yr)
- `surface_dx`: ice surface gradient in x (∂s/∂x, unitless)
- `surface_dy`: ice surface gradient in y (∂s/∂y, unitless)
"""
struct ElevationState{AA,ACX,ACY}
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
struct TopographicState{B,M,ACX,ACY}
    mask::TopographyMasks{B}
    distance::DistanceState{M}
    fraction::FractionState{M}
    thickness::ThicknessState{M}
    massbalance::MassBalanceState{M}
    elevation::ElevationState{M,ACX,ACY}
end
Adapt.@adapt_structure TopographicState

function TopographicState(grid::RegularGrid)
    backend = KernelAbstractions.get_backend(grid.x)
    T = eltype(grid.x)
    (; nx, ny) = grid
    b = KernelAbstractions.zeros(backend, Bool, nx, ny)
    m = KernelAbstractions.zeros(backend, T, nx, ny)
    return TopographicState(
        TopographyMasks([copy(b) for _ = 1:7]...),
        DistanceState([copy(m) for _ = 1:2]...),
        FractionState([copy(m) for _ = 1:1]...),
        ThicknessState([copy(m) for _ = 1:6]...),
        MassBalanceState([copy(m) for _ = 1:10]...),
        ElevationState([copy(m) for _ = 1:9]...),
    )
end

"""
$(TYPEDSIGNATURES)

Build a [`TopographicState`](@ref) of Chmy `Field`s on `grid`. All fields are
depth-integrated, so they live on `grid.grid2d` at the `aa` node — except the surface
gradients, which live on the `acx`/`acy` faces. Element type follows the grid; `halo` is the
ghost-cell width on `x`/`y` — every field here is on `grid.grid2d`, so `z` gets none
(see [`_halo2d`](@ref)).

`diagnostics = false` (the default) allocates the output-only fields on a one-cell grid
instead of the full one — see the "Diagnostic (output-only) fields" note above. In this
struct that covers both `DistanceState` fields, nine of the ten `MassBalanceState` terms
(all but `net`), `ThicknessState`'s `ice_ref`/`ice_grounded`/`sediment`, and
`ElevationState`'s `bed_ref`/`bed_stddev`/`surface_dt`: 17 of 35 fields.
"""
function TopographicState(grid::StaggeredGrid; halo = 1, diagnostics = false)
    (; arch) = grid
    g = grid.grid2d
    T = eltype(g)
    d = _degenerate_grid(grid)
    h = _halo2d(halo)
    b() = _field(arch, g, NODE_AA, Bool, h)
    aa() = _field(arch, g, NODE_AA, T, h)
    acx() = _field(arch, g, NODE_ACX, T, h)
    acy() = _field(arch, g, NODE_ACY, T, h)
    # @diagnostic — Necessary ✗ / Used ✗ in `docs/src/variables.md`
    aa_d() = _maybe_field(arch, g, d, NODE_AA, T, h, diagnostics)
    return TopographicState(
        TopographyMasks(ntuple(_ -> b(), 7)...),
        DistanceState(aa_d(), aa_d()),                       # distance_to_margin/grline
        FractionState(aa()),
        # ice, ice_ref, ice_dt, ice_effective, ice_grounded, sediment
        ThicknessState(aa(), aa_d(), aa(), aa(), aa_d(), aa_d()),
        # base, base_floating, base_grounded, calving_floating, calving_grounded,
        # discharge, front, net, surface, surface_ref — only `net` is live
        MassBalanceState(ntuple(_ -> aa_d(), 7)..., aa(), aa_d(), aa_d()),
        # base, bed, bed_ref, bed_stddev, seasurface, surface, surface_dt, surface_dx/dy
        ElevationState(aa(), aa(), aa_d(), aa_d(), aa(), aa(), aa_d(), acx(), acy()),
    )
end

###############################################################

"""
$(TYPEDSIGNATURES)

## Fields
- `surface`: ice surface velocity (m/yr)
- `thickness`: ice thickness (m)
"""
struct MechanicTopographyState{AA}
    surface::AA
    thickness::AA
end
Adapt.@adapt_structure MechanicTopographyState

"""
$(TYPEDSIGNATURES)

Material properties feeding the momentum balance. `viscosity`/`viscosity_depthaveraged` are
the 3D and depth-averaged ice viscosity. `rate_factor`/`rate_factor_depthaveraged` are
prescribed inputs, not derived from temperature (`ThermodynamicState` is a separate,
unconnected sibling; see `pagos-roadmap/PT-autotune.md`, Phase 1); they exist only for the
pseudo-transient solver's Glen-law viscosity continuation
([`GlenViscosityContinuation`](@ref)), and `rate_factor` is filled by broadcasting
`rate_factor_depthaveraged` down the column in the isothermal case. `viscosity_integral_1`/
`viscosity_integral_2` are DIVA's generalized viscosity integrals `F_m = ∫_b^s
(1/µ)((s-z)/H)^m dz` (Robinson et al. 2022, Eq. 15), written by `viscosity_integrals!` and
kept here next to the viscosity they integrate rather than in `FrictionState`; zero on any
state that never runs DIVA.

## Fields
- `viscosity_depthaveraged`: depth-averaged ice viscosity (Pa yr)
- `viscosity`: 3D ice viscosity (Pa yr)
- `rate_factor_depthaveraged`: depth-averaged ice rate factor (Pa⁻¹ yr⁻¹)
- `rate_factor`: 3D ice rate factor (Pa⁻¹ yr⁻¹)
- `viscosity_integral_1`: first generalized viscosity integral
- `viscosity_integral_2`: second generalized viscosity integral
"""
struct MechanicMaterialState{AA2,AA3}
    viscosity_depthaveraged::AA2
    viscosity::AA3
    rate_factor_depthaveraged::AA2
    rate_factor::AA3
    viscosity_integral_1::AA2
    viscosity_integral_2::AA2
end
Adapt.@adapt_structure MechanicMaterialState

# beta, beta_eff, c_bed live at `aa`, where the effective pressure and basal velocity
# magnitude they depend on are evaluated. Whether β on the velocity faces (β_acx/acy)
# becomes stored fields or an inline `lerp` in the momentum kernel is a Phase 3 decision
# (see `pagos-roadmap/chmy.md`).
"""
$(TYPEDSIGNATURES)

## Fields
- `beta`: basal friction coefficient (Pa yr m⁻¹)
- `beta_eff`: effective basal friction coefficient (Pa yr m⁻¹)
- `c_bed`: basal yield stress (Pa)
"""
struct FrictionState{AA}
    beta::AA
    beta_eff::AA
    c_bed::AA
end
Adapt.@adapt_structure FrictionState

# The flux q = Hū must sit on the cell faces the continuity equation's divergence reads
# (acx/acy) for `∂H/∂t = -divg(q)` to stay conservative — a cell-centred `aa` flux would
# need re-staggering inside the divergence. `grline` is a grounding-line diagnostic, not a
# continuity term, so it stays at `aa`.
"""
$(TYPEDSIGNATURES)

## Fields
- `x`: ice flux in x
- `y`: ice flux in y
- `grline`: grounding line diagnostic
"""
struct FluxState{ACX2,ACY2,AA2}
    x::ACX2
    y::ACY2
    grline::AA2
end
Adapt.@adapt_structure FluxState

"""
$(TYPEDSIGNATURES)

## Fields
- `driving_x`: driving stress in x (Pa)
- `driving_y`: driving stress in y (Pa)
- `base_x`: basal stress in x (Pa)
- `base_y`: basal stress in y (Pa)
- `base_vertical`: basal stress in z (Pa)
- `membrane_xx`: depth-integrated membrane stress in xx (Pa m)
- `membrane_xy`: depth-integrated membrane stress in xy (Pa m)
- `membrane_yy`: depth-integrated membrane stress in yy (Pa m)
- `xx`: stress tensor component σ_xx (Pa)
- `xy`: stress tensor component σ_xy (Pa)
- `xz`: stress tensor component σ_xz (Pa)
- `yy`: stress tensor component σ_yy (Pa)
- `yz`: stress tensor component σ_yz (Pa)
- `zz`: stress tensor component σ_zz (Pa)
- `effective`: effective stress (Pa)
- `lateral`: lateral stress (Pa)
- `eigenvalue_1`: first principal stress (Pa)
- `eigenvalue_2`: second principal stress (Pa)

`yx`/`zx`/`zy` are not stored: the tensor is symmetric, so they are exposed as
non-allocating aliases of `xy`/`xz`/`yz` via `getproperty` (see below `StressState`).
"""
struct StressState{ACX2,ACY2,AA2,AB2,AA3,AB3,ACXZ3,ACYZ3}
    driving_x::ACX2
    driving_y::ACY2
    base_x::ACX2
    base_y::ACY2
    base_vertical::AA2

    # Depth-integrated membrane stress the SSA/DIVA momentum balance differentiates:
    # membrane_xx = 2µ̄H(2ūx + v̄y), membrane_xy = µ̄H(ūy + v̄x),
    # membrane_yy = 2µ̄H(ūx + 2v̄y) (Robinson et al. 2022, Eq. 14) — `aa` for the normal
    # components, `ab` for the shear one, where the gradients that build them live.
    membrane_xx::AA2
    membrane_xy::AB2
    membrane_yy::AA2

    # Symmetric: `yx`/`zx`/`zy` are not fields here, but `getproperty` below aliases them
    # to `xy`/`xz`/`yz` so consumers can still spell out the full tensor.
    xx::AA3
    xy::AB3
    xz::ACXZ3
    yy::AA3
    yz::ACYZ3
    zz::AA3
    effective::AA3
    lateral::AA3
    eigenvalue_1::AA3
    eigenvalue_2::AA3
end
Adapt.@adapt_structure StressState

function Base.getproperty(s::StressState, name::Symbol)
    name === :yx && return getfield(s, :xy)
    name === :zx && return getfield(s, :xz)
    name === :zy && return getfield(s, :yz)
    return getfield(s, name)
end

"""
$(TYPEDSIGNATURES)

## Fields
- `xx`: strain rate tensor component ε̇_xx (yr⁻¹)
- `xy`: strain rate tensor component ε̇_xy (yr⁻¹)
- `xz`: strain rate tensor component ε̇_xz (yr⁻¹)
- `yy`: strain rate tensor component ε̇_yy (yr⁻¹)
- `yz`: strain rate tensor component ε̇_yz (yr⁻¹)
- `zz`: strain rate tensor component ε̇_zz (yr⁻¹)
- `effective`: effective strain rate (yr⁻¹)
- `effective_depthaveraged`: depth-averaged effective strain rate (yr⁻¹)

`yx`/`zx`/`zy` are not stored: the tensor is symmetric, so they are exposed as
non-allocating aliases of `xy`/`xz`/`yz` via `getproperty` (see below `StrainRateState`).
"""
struct StrainRateState{AA2,AA3,AB3,ACXZ3,ACYZ3}
    # Symmetric: `yx`/`zx`/`zy` are not fields here, but `getproperty` below aliases them
    # to `xy`/`xz`/`yz` so consumers can still spell out the full tensor.
    xx::AA3
    xy::AB3
    xz::ACXZ3
    yy::AA3
    yz::ACYZ3
    zz::AA3

    # `effective` and `effective_depthaveraged` are different quantities, not one value at
    # two resolutions. `effective` is DIVA's (Robinson et al. 2022, Eq. 13): it carries the
    # vertical-shear terms ¼(u_z² + v_z²), differs layer by layer, and feeds the 3D
    # `material.viscosity`. `effective_depthaveraged` is the SSA one (Eq. 12), the same
    # expression with the shear terms dropped, feeding `material.viscosity_depthaveraged`.
    effective::AA3
    effective_depthaveraged::AA2
end
Adapt.@adapt_structure StrainRateState

function Base.getproperty(s::StrainRateState, name::Symbol)
    name === :yx && return getfield(s, :xy)
    name === :zx && return getfield(s, :xz)
    name === :zy && return getfield(s, :yz)
    return getfield(s, name)
end

struct VelocityState{ACX2,ACY2,AA2,AB2,ACX3,ACY3,AA3,AB3,AAZ3,ACXZ3,ACYZ3}
    depthaverage_x::ACX2
    depthaverage_y::ACY2
    depthaverage_x_dx::AA2
    depthaverage_x_dy::AB2
    depthaverage_x_dz::ACX2
    depthaverage_y_dx::AB2
    depthaverage_y_dy::AA2
    depthaverage_y_dz::ACY2

    base_x::ACX2
    base_y::ACY2
    base_norm::AA2

    surface_x::ACX2
    surface_y::ACY2
    surface_norm::AA2

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
struct MechanicState{ACX2,ACY2,AA2,AB2,ACX3,ACY3,AA3,AB3,AAZ3,ACXZ3,ACYZ3}
    friction::FrictionState{AA2}
    flux::FluxState{ACX2,ACY2,AA2}

    topography::MechanicTopographyState{AA2}
    material::MechanicMaterialState{AA2,AA3}
    strainrate::StrainRateState{AA2,AA3,AB3,ACXZ3,ACYZ3}
    stress::StressState{ACX2,ACY2,AA2,AB2,AA3,AB3,ACXZ3,ACYZ3}
    velocity::VelocityState{ACX2,ACY2,AA2,AB2,ACX3,ACY3,AA3,AB3,AAZ3,ACXZ3,ACYZ3}
end
Adapt.@adapt_structure MechanicState

function MechanicState(grid::RegularGrid)
    backend = KernelAbstractions.get_backend(grid.x)
    T = eltype(grid.x)
    (; nx, ny, nz) = grid
    m2() = KernelAbstractions.zeros(backend, T, nx, ny)
    m3() = KernelAbstractions.zeros(backend, T, nx, ny, nz)
    return MechanicState(
        FrictionState(m2(), m2(), m2()),                    # beta, beta_eff, c_bed
        FluxState(m2(), m2(), m2()),                       # flux x, y, grline
        MechanicTopographyState(m2(), m2()),               # surface, thickness
        # viscosity_depthaveraged, viscosity, rate_factor_depthaveraged, rate_factor, F₁, F₂
        MechanicMaterialState(m2(), m3(), m2(), m3(), m2(), m2()),
        StrainRateState(ntuple(_ -> m3(), 7)..., m2()),   # 7 column + effective_depthaveraged
        StressState(ntuple(_ -> m2(), 8)..., ntuple(_ -> m3(), 10)...),  # 8 depth-integrated + 10 column
        VelocityState(ntuple(_ -> m2(), 14)..., ntuple(_ -> m3(), 13)...),  # 14 depth-averaged + 13 column
    )
end

"""
$(TYPEDSIGNATURES)

Build a [`MechanicState`](@ref) of Chmy `Field`s on `grid`. Depth-integrated fields are
built on `grid.grid2d`, column fields on `grid.grid`; each is placed at the node class
its physics dictates (see the table at the top of `src/api/state.jl`). Element type
follows the grid; `halo` is the ghost-cell width on `x`/`y`/(`z` for column fields) — most
depth-integrated fields get no `z` ghost (see [`_halo2d`](@ref)).

!!! note "The `ACX2`/`ACY2` fields are the exception, and keep a full `z` ghost"
    `flux.x`/`y`, `stress.driving_x`/`y`, `stress.base_x`/`y`, and eight fields in
    [`VelocityState`](@ref) (`depthaverage_x`/`y`, `depthaverage_x_dz`/`y_dz`, `base_x`/`y`,
    `surface_x`/`y`) all carry `MechanicState`'s own `ACX2`/`ACY2` type parameter — the
    *outer* struct's, reused across `flux`/`stress`/`velocity` rather than each owning its
    own — so they share one halo whether or not a given field needs the full one.
    `velocity.depthaverage_x`/`y` are `bc!`'d every PT iteration
    (`pseudo_transient!`), and Chmy's boundary-condition kernel sweeps every axis *other*
    than the one being filled at the grid's own `size(·, Vertex()) .+ 2` — so even filling
    the `x`/`y` boundary touches a `z` "ghost" that a `(h, h, 0)`-halo field does not have
    (a `BoundsError` under `--check-bounds=yes`; silent, intermittent corruption without
    it — see `pagos-roadmap/chmy-issues.md`). These 14 fields are therefore exempt.

`diagnostics = false` (the default) allocates the output-only fields on a one-cell grid
instead of the full one — see the "Diagnostic (output-only) fields" note above. Here that is
`flux.grline`, `stress.base_vertical`, `stress.lateral`, `stress.eigenvalue_1/2` and
`velocity.norm`. The four column fields among them are the expensive ones: at 1522×1522 with
`nz = 11` they are ~140 MB each.

Note this does **not** cover the full Cauchy stress tensor (`stress.xx`…`zz`) or
`stress.effective`: those are diagnostic in the sense of never being read back into the
solve, but `src/mechanics/stress.jl` does actively *write* them, so making them degenerate
would break that path. Gating those on an opt-in is worth a further ~1.4 GiB at 4 km and is
left as a separate change.
"""
function MechanicState(grid::StaggeredGrid; halo = 1, diagnostics = false)
    (; arch) = grid
    g2, g3 = grid.grid2d, grid.grid
    T = eltype(g3)
    d = _degenerate_grid(grid)
    h2 = _halo2d(halo)
    # @diagnostic — Necessary ✗ / Used ✗ in `docs/src/variables.md`
    aa2_d() = _maybe_field(arch, g2, d, NODE_AA, T, h2, diagnostics)
    aa3_d() = _maybe_field(arch, g3, d, NODE_AA, T, halo, diagnostics)
    aa2() = _field(arch, g2, NODE_AA, T, h2)
    ab2() = _field(arch, g2, NODE_AB, T, h2)
    # NOT `h2`: `MechanicState`'s own type parameters unify `ACX2`/`ACY2` across `flux`,
    # `stress` and `velocity` (they're the *outer* struct's params, reused — not each
    # sub-state's own), so every field built from `acx2()`/`acy2()` shares one halo whether
    # or not that particular field needs it. `velocity.depthaverage_x`/`y` are `bc!`'d every
    # PT iteration (`pseudo_transient!`), and Chmy's boundary-condition kernel sweeps every
    # *other* axis at `size(grid, Vertex()) .+ 2` regardless of the target field's own halo
    # — so filling even the `x`/`y` boundary reads/writes a `z` "ghost" a `(h,h,0)` field
    # does not have. Confirmed with `--check-bounds=yes`: `BoundsError: ... at index [2, 3,
    # 0]`, non-deterministically past `@inbounds` on a normal build. See
    # `pagos-roadmap/chmy-issues.md` for the upstream report. Costs 14 of `MechanicState`'s
    # ~35 depth-integrated fields the full saving (flux.x/y, stress.driving_x/y,
    # stress.base_x/y, and 8 in `VelocityState` — see its docstring); the `aa2`/`ab2` fields
    # above, and everything in `TopographicState`/`ThermodynamicState`/`MaterialState`
    # (never `bc!`'d), are unaffected.
    acx2() = _field(arch, g2, NODE_ACX, T, halo)
    acy2() = _field(arch, g2, NODE_ACY, T, halo)
    aa3() = _field(arch, g3, NODE_AA, T, halo)
    acx3() = _field(arch, g3, NODE_ACX, T, halo)
    acy3() = _field(arch, g3, NODE_ACY, T, halo)
    ab3() = _field(arch, g3, NODE_AB, T, halo)
    aaz3() = _field(arch, g3, NODE_AA_AC, T, halo)
    acxz3() = _field(arch, g3, NODE_ACX_AC, T, halo)
    acyz3() = _field(arch, g3, NODE_ACY_AC, T, halo)
    return MechanicState(
        FrictionState(aa2(), aa2(), aa2()),
        FluxState(acx2(), acy2(), aa2_d()),                # flux x, y, grline(@diagnostic)
        MechanicTopographyState(aa2(), aa2()),
        # viscosity_depthaveraged, viscosity, rate_factor_depthaveraged, rate_factor, F₁, F₂
        MechanicMaterialState(aa2(), aa3(), aa2(), aa3(), aa2(), aa2()),
        StrainRateState(
            aa3(),
            ab3(),
            acxz3(),
            aa3(),         # xx, xy, xz, yy
            acyz3(),
            aa3(),            # yz, zz
            aa3(),
            aa2(),                                # effective, effective_depthaveraged
        ),
        StressState(
            acx2(),
            acy2(),
            acx2(),
            acy2(),
            aa2_d(),     # driving_x/y, base_x/y, base_vertical(@diagnostic)
            aa2(),
            ab2(),
            aa2(),                         # membrane_xx, membrane_xy, membrane_yy
            aa3(),
            ab3(),
            acxz3(),
            aa3(),         # xx, xy, xz, yy
            acyz3(),
            aa3(),            # yz, zz
            aa3(),            # effective
            aa3_d(),
            aa3_d(),
            aa3_d(),          # lateral, eigenvalue_1/2 (@diagnostic)
        ),
        VelocityState(
            acx2(),
            acy2(),
            aa2(),
            ab2(),
            acx2(),        # depthaverage_x, depthaverage_y, depthaverage_x_dx/dy/dz
            ab2(),
            aa2(),
            acy2(),                        # depthaverage_y_dx/dy/dz
            acx2(),
            acy2(),
            aa2(),                       # base_x, base_y, base_norm
            acx2(),
            acy2(),
            aa2(),                       # surface_x, surface_y, surface_norm
            acx3(),
            acy3(),
            aaz3(),                      # x, y, z
            aa3(),
            ab3(),
            acxz3(),                       # x_dx, x_dy, x_dz
            ab3(),
            aa3(),
            acyz3(),                       # y_dx, y_dy, y_dz
            acxz3(),
            acyz3(),
            aa3(),              # z_dx, z_dy, z_dz
            aa3_d(),            # norm (@diagnostic)
        ),
    )
end

###############################################################

struct TemperatureState{AA3,AA2}
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
struct ThermodynamicState{AA3,AA2}
    temperature::TemperatureState{AA3,AA2}
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
    T = eltype(grid.x)
    (; nx, ny, nz) = grid
    m() =
        nz > 1 ? KernelAbstractions.zeros(backend, T, nx, ny, nz) :
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
grid; `halo` is the ghost-cell width on `x`/`y`/(`z` for column fields) — the
depth-integrated fields get no `z` ghost (see [`_halo2d`](@ref)).
"""
function ThermodynamicState(grid::StaggeredGrid; halo = 1, diagnostics = false)
    (; arch) = grid
    g2, g3 = grid.grid2d, grid.grid
    T = eltype(g3)
    d = _degenerate_grid(grid)
    h2 = _halo2d(halo)
    aa2() = _field(arch, g2, NODE_AA, T, h2)
    aa3() = _field(arch, g3, NODE_AA, T, halo)
    # @diagnostic — Necessary ✗ / Used ✗ in `docs/src/variables.md`. Every *other* field here
    # is Used ✗ too (the heat equation is not implemented yet) but Necessary ✓, so it stays
    # full size: it is waiting on a solver, not on a decision.
    aa2_d() = _maybe_field(arch, g2, d, NODE_AA, T, h2, diagnostics)
    aa3_d() = _maybe_field(arch, g3, d, NODE_AA, T, halo, diagnostics)
    return ThermodynamicState(
        TemperatureState(aa3(), aa2(), aa3(), aa3(), aa3()),
        EnthalpyState(aa3(), aa3()),
        aa3(),
        aa3(),
        aa3_d(),                     # ice_water_content, heat_strain_internal, _dt(@diagnostic)
        aa2(),
        aa2(),
        aa2(),
        aa2(),                # heat_base_friction, heatflux_base_ice/bedrock/geothermal
        aa2_d(),
        aa2_d(),
        aa2_d(),              # thickness_waterlayer(_dt), coldtemperate_interface (@diagnostic)
        aa3(),
        aa3(),                              # specific_heat_capacity_ice, heat_conductivity_ice
    )
end

###############################################################

"""
$(TYPEDSIGNATURES)

State variables for the material component, all at the `aa` node: `AA2` are the
depth-averaged and depth-integrated viscosities, `AA3` the viscosity of the ice column.
"""
struct MaterialState{AA2,AA3}
    eta_depth_averaged::AA2
    eta_depth_integrated::AA2
    eta_ice::AA3
end
Adapt.@adapt_structure MaterialState

function MaterialState(grid::RegularGrid)
    backend = KernelAbstractions.get_backend(grid.x)
    T = eltype(grid.x)
    (; nx, ny, nz) = grid
    m2 = KernelAbstractions.zeros(backend, T, nx, ny)
    m3 = KernelAbstractions.zeros(backend, T, nx, ny, nz)
    return MaterialState(copy(m2), copy(m2), m3)
end

"""
$(TYPEDSIGNATURES)

Build a [`MaterialState`](@ref) of Chmy `Field`s on `grid`: the depth-averaged and
depth-integrated viscosities on `grid.grid2d`, the column viscosity on `grid.grid`, all
at the `aa` node. Element type follows the grid; `halo` is the ghost-cell width on
`x`/`y`/(`z` for `eta_ice`) — the two depth-integrated fields get no `z` ghost (see
[`_halo2d`](@ref)).

`diagnostics = false` (the default) allocates `eta_depth_integrated` on a one-cell grid —
see the "Diagnostic (output-only) fields" note above. It is the depth-*integral*
`∫ η dz`, distinct from the depth-*average* `η̄` that is live, and nothing reads it.
"""
function MaterialState(grid::StaggeredGrid; halo = 1, diagnostics = false)
    (; arch) = grid
    T = eltype(grid.grid)
    h2 = _halo2d(halo)
    aa2() = _field(arch, grid.grid2d, NODE_AA, T, h2)
    # @diagnostic — Necessary ✗ / Used ✗ in `docs/src/variables.md`
    aa2_d() = _maybe_field(arch, grid.grid2d, _degenerate_grid(grid), NODE_AA, T, h2, diagnostics)
    return MaterialState(aa2(), aa2_d(), _field(arch, grid.grid, NODE_AA, T, halo))
end

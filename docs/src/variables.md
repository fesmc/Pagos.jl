# [Variables](@id variables)

This page inventories the fields carried by Pagos' state structs (`src/api/state.jl`), grouped by the top-level model component that owns them: [Topography](@ref topography), [Mechanics](@ref) (dynamics), [Thermodynamics](@ref thermodynamics_state), and [Material](@ref material). For the grid location (`aa`/`acx`/`acy`/`ab`, and the `2`/`3` depth-integrated/column suffix) each field is built on, see [Staggered grids](@ref staggered_grids).

**Columns**

- **Variable** — short English name.
- **Symbol** — the math notation used elsewhere in these docs.
- **Dimension** — `2D` for a depth-integrated/horizontal-only field; `3D` for a field that always carries a vertical layer index; `2D/3D` for a column field whose vertical extent depends on a modelling choice — it collapses to a single layer (`nz == 1`) under a depth-averaged momentum balance (SIA/SSA) or a single-layer thermodynamics setup, and has `nz > 1` layers otherwise (see [the `2`/`3` suffix](@ref node_type_parameters)).
- **Instances** — which state struct(s) carry this field. The same physical quantity sometimes lives in more than one struct because different components don't necessarily share a grid — e.g. viscosity is stored once in `MaterialState` (the material component's own grid) and again in `MechanicMaterialState` (the copy the momentum solver actually reads).
- **Necessary** — ✓ if the variable feeds the prognostic update of ice thickness, velocity, temperature, or viscosity (directly, or as an intermediate on that path); ✗ if it is a diagnostic computed for output/analysis and never read back into the next step. This reflects the variable's intended role in the model design, independent of whether that role is already implemented.
- **Used** — ✓ if some function in `src/` currently reads and/or writes this field; ✗ if the field is allocated in the state struct but nothing in `src/` touches it yet. A field can be necessary but currently unused (an as-designed dependency on a component that isn't wired up yet, e.g. most of `ThermodynamicState`), and a field can be used without being internally computed — several fields (e.g. `mask.is_grounded`, `mask.is_ice_allowed`, `massbalance.net`) are deliberately caller-supplied forcing/inputs that live code reads, rather than something a solver derives.

## Topography

| Variable | Symbol | Dimension | Instances | Necessary | Used |
|---|---|---|---|---|---|
| ice mask | $\mathbb{1}_{\mathrm{ice}}$ | 2D | `TopographyMasks` | ✓ | ✓ |
| ice-neighbour mask | $\mathbb{1}_{\mathrm{nbr}}$ | 2D | `TopographyMasks` | ✓ | ✓ |
| allowed-ice mask | $\mathbb{1}_{\mathrm{allowed}}$ | 2D | `TopographyMasks` | ✓ | ✓ |
| grounded mask | $\mathbb{1}_{\mathrm{grounded}}$ | 2D | `TopographyMasks` | ✓ | ✓ |
| floating mask | $\mathbb{1}_{\mathrm{floating}}$ | 2D | `TopographyMasks` | ✓ | ✗ |
| margin mask | $\mathbb{1}_{\mathrm{margin}}$ | 2D | `TopographyMasks` | ✓ | ✗ |
| momentum-solved mask | $\mathbb{1}_{\mathrm{solved}}$ | 2D | `TopographyMasks` | ✓ | ✓ |
| distance to margin | $d_{\mathrm{margin}}$ | 2D | `DistanceState` | ✗ | ✗ |
| distance to grounding line | $d_{\mathrm{gl}}$ | 2D | `DistanceState` | ✗ | ✗ |
| grounded fraction | $f_{\mathrm{g}}$ | 2D | `FractionState` | ✓ | ✗ |
| net mass balance | $\dot m_{\mathrm{net}}$ | 2D | `MassBalanceState` | ✓ | ✓ |
| basal mass balance | $\dot m_{\mathrm{base}}$ | 2D | `MassBalanceState` | ✗ | ✗ |
| sub-shelf mass balance | $\dot m_{\mathrm{base,fl}}$ | 2D | `MassBalanceState` | ✗ | ✗ |
| grounded-base mass balance | $\dot m_{\mathrm{base,gr}}$ | 2D | `MassBalanceState` | ✗ | ✗ |
| floating calving loss | $\dot m_{\mathrm{calv,fl}}$ | 2D | `MassBalanceState` | ✗ | ✗ |
| grounded calving loss | $\dot m_{\mathrm{calv,gr}}$ | 2D | `MassBalanceState` | ✗ | ✗ |
| grounding-line discharge | $\dot m_{\mathrm{disch}}$ | 2D | `MassBalanceState` | ✗ | ✗ |
| frontal mass balance | $\dot m_{\mathrm{front}}$ | 2D | `MassBalanceState` | ✗ | ✗ |
| surface mass balance | $\dot m_{\mathrm{srf}}$ | 2D | `MassBalanceState` | ✗ | ✗ |
| reference surface mass balance | $\dot m_{\mathrm{srf,ref}}$ | 2D | `MassBalanceState` | ✗ | ✗ |
| ice thickness | $H$ | 2D | `ThicknessState`, `MechanicTopographyState` | ✓ | ✓ |
| thickness tendency | $\partial H/\partial t$ | 2D | `ThicknessState` | ✓ | ✓ |
| effective ice thickness | $H_{\mathrm{eff}}$ | 2D | `ThicknessState` | ✓ | ✗ |
| reference ice thickness | $H_{\mathrm{ref}}$ | 2D | `ThicknessState` | ✗ | ✗ |
| grounded ice thickness | $H_{\mathrm{grounded}}$ | 2D | `ThicknessState` | ✗ | ✗ |
| sediment thickness | $H_{\mathrm{sed}}$ | 2D | `ThicknessState` | ✗ | ✗ |
| bed elevation | $z_{\mathrm{b}}$ | 2D | `ElevationState` | ✓ | ✗ |
| sea level | $z_{\mathrm{sl}}$ | 2D | `ElevationState` | ✓ | ✗ |
| surface elevation | $z_{\mathrm{srf}}$ | 2D | `ElevationState`, `MechanicTopographyState` | ✓ | ✓ |
| ice base elevation | $z_{\mathrm{base}}$ | 2D | `ElevationState` | ✓ | ✗ |
| surface slope (x) | $\partial z_{\mathrm{srf}}/\partial x$ | 2D | `ElevationState` | ✓ | ✓ |
| surface slope (y) | $\partial z_{\mathrm{srf}}/\partial y$ | 2D | `ElevationState` | ✓ | ✓ |
| reference bed elevation | $z_{\mathrm{b,ref}}$ | 2D | `ElevationState` | ✗ | ✗ |
| bed roughness | $\sigma_{z_{\mathrm{b}}}$ | 2D | `ElevationState` | ✗ | ✗ |
| surface elevation tendency | $\partial z_{\mathrm{srf}}/\partial t$ | 2D | `ElevationState` | ✗ | ✗ |

Every topography field is depth-integrated (built on `grid2d`), so this component has no `2D/3D` rows.

All seven masks *gate* computations by design and are necessary in that sense, but only four are currently live: `is_ice`/`is_ice_neighbour` are computed by `icemasks!`, and `is_momentum_solved` by `momentum_mask!`; `is_ice_allowed` and `is_grounded` are read (not computed — per the comment at `src/topography/masks.jl:276`, both are deliberately caller-supplied, like an imposed domain boundary or an externally determined grounding state). `is_floating` and `is_margin` are allocated but not referenced anywhere else in `src/` yet. `distance_to_margin`/`distance_to_grline` are likewise unused placeholders. `fraction_grounded` is necessary once wired (Coulomb-type friction laws take a grounded-fraction argument, `src/mechanics/basal_friction.jl`), but the field itself isn't currently populated or read anywhere — `src/topography/subgrid.jl`'s grounded-fraction functions exist but aren't yet wired to write into it.

Only `net` is currently read by the continuity solver (`src/topography/advection.jl`) — it is the caller-supplied forcing that $\partial H/\partial t = -\nabla\cdot q + \dot m_{\mathrm{net}}$ advances thickness with (see `test/api/advection.jl`). The disaggregated per-process terms (basal, calving, discharge, frontal, surface) are allocated as a place to store a component-wise mass-balance breakdown for output, but nothing in `src/` currently sums them into `net`, so treat them as diagnostic *and* unused for now.

$H$ is the core prognostic variable of the topography component, read by `icemasks!` and advanced by `advect!`; `MechanicTopographyState.thickness` is the copy the momentum solver reads on its own grid, kept in sync separately. $H_{\mathrm{eff}} = H/f_{\mathrm{ice}}$ is necessary by design (`src/topography/height.jl` defines it, feeding surface elevation, height-above-floatation, and calving-rate parameterizations), but — like `fraction_grounded` above — the dedicated state field is not yet populated by any of those functions. `ice_ref`, `ice_grounded`, and `sediment` are unused placeholders.

$z_{\mathrm{b}}$, $z_{\mathrm{sl}}$ and $z_{\mathrm{base}}$ are necessary to set the flotation criterion the grounded/floating masks are built from, but — again — none of the three is currently referenced as a state field outside `src/api/state.jl`; the flotation logic in `src/topography/height.jl` exists only as free functions over plain arguments, not yet wired to these fields. `surface`, by contrast, is live: `surface_gradient!` reads it and writes the surface slope, which directly builds the driving stress $\tau_{\mathrm{d}}$. `bed_ref` and `bed_stddev` are unused placeholders (the latter reserved for a future bed-roughness-dependent friction law); `surface_dt` is allocated but unused.

## Mechanics

| Variable | Symbol | Dimension | Instances | Necessary | Used |
|---|---|---|---|---|---|
| basal friction coefficient | $\beta$ | 2D | `FrictionState` | ✓ | ✓ |
| effective friction coefficient | $\beta_{\mathrm{eff}}$ | 2D | `FrictionState` | ✓ | ✓ |
| bed strength | $c_{\mathrm{b}}$ | 2D | `FrictionState` | ✓ | ✗ |
| ice flux | $q_x, q_y$ | 2D | `FluxState` | ✓ | ✓ |
| grounding-line flux | $q_{\mathrm{gl}}$ | 2D | `FluxState` | ✗ | ✗ |
| driving stress | $\tau_{\mathrm{d},x}, \tau_{\mathrm{d},y}$ | 2D | `StressState` | ✓ | ✓ |
| basal stress | $\tau_{\mathrm{b},x}, \tau_{\mathrm{b},y}$ | 2D | `StressState` | ✓ | ✓ |
| membrane stress | $\tau_{\mathrm{m},xx}, \tau_{\mathrm{m},xy}, \tau_{\mathrm{m},yy}$ | 2D | `StressState` | ✓ | ✓ |
| Cauchy stress tensor | $\sigma_{ij}$ | 2D/3D | `StressState` | ✗ | ✓ |
| effective stress | $\sigma_{\mathrm{e}}$ | 2D/3D | `StressState` | ✗ | ✓ |
| vertical basal stress | $\tau_{\mathrm{b},z}$ | 2D | `StressState` | ✗ | ✗ |
| lateral stress | $\sigma_{\mathrm{lat}}$ | 2D/3D | `StressState` | ✗ | ✗ |
| principal stresses | $\sigma_1, \sigma_2$ | 2D/3D | `StressState` | ✗ | ✗ |
| strain-rate tensor | $\dot\varepsilon_{ij}$ | 2D/3D | `StrainRateState` | ✓ | ✓ |
| 3D effective strain rate | $\dot\varepsilon_{\mathrm{e}}(z)$ | 2D/3D | `StrainRateState` | ✓ | ✓ |
| depth-averaged effective strain rate | $\bar{\dot\varepsilon}_{\mathrm{e}}$ | 2D | `StrainRateState` | ✓ | ✓ |
| depth-averaged velocity | $\bar u, \bar v$ | 2D | `VelocityState` | ✓ | ✓ |
| depth-averaged velocity gradient | $\partial \bar u_i/\partial x_j$ | 2D | `VelocityState` | ✓ | ✓ |
| basal velocity | $u_{\mathrm{b}}, v_{\mathrm{b}}, \lvert v_{\mathrm{b}}\rvert$ | 2D | `VelocityState` | ✓ | ✓ |
| 3D velocity | $u(z), v(z), w(z)$ | 2D/3D | `VelocityState` | ✓ | ✓ |
| velocity-gradient tensor | $\partial u_i/\partial x_j$ | 2D/3D | `VelocityState` | ✓ | ✓ |
| surface velocity | $u_{\mathrm{s}}, v_{\mathrm{s}}, \lvert v_{\mathrm{s}}\rvert$ | 2D | `VelocityState` | ✗ | ✓ |
| 3D velocity magnitude | $\lvert v(z)\rvert$ | 2D/3D | `VelocityState` | ✗ | ✗ |
| DIVA viscosity integrals | $F_1, F_2$ | 2D | `MechanicMaterialState` | ✓ | ✓ |
| column rate factor | $A(z)$ | 2D/3D | `MechanicMaterialState` | ✓ | ✓ |
| depth-averaged rate factor | $\bar A$ | 2D | `MechanicMaterialState` | ✓ | ✓ |

`MechanicTopographyState` and viscosity ($\eta(z)$, $\bar\eta$) are documented in [Topography](@ref topography) and [Material](@ref material) respectively, since they are the same physical quantities as their counterparts there, just a per-component copy on the mechanics grid.

All three friction fields are necessary — they enter the basal-stress term of the momentum balance, and $\beta_{\mathrm{eff}}$ is the regularized form the linear solve actually consumes — but `c_bed` is not yet referenced anywhere in `src/` (reserved for a Coulomb-type friction law that isn't wired up yet). $q = H\bar u$ is the term differenced in $\partial H/\partial t = -\nabla\cdot q$; `grline` is explicitly a grounding-line diagnostic, not a continuity term, and is currently unused (see the comment at `src/api/state.jl:249`).

$\tau_{\mathrm{d}}$ is the driving-stress forcing and $\tau_{\mathrm{m}}$ the depth-integrated membrane stress the SSA/DIVA momentum balance actually differentiates ([robinson_comparison_2022](@citet) Eq. 14) — both are live. Viscosity in this architecture is built from effective *strain rate*, not effective stress, so the full 9-component Cauchy tensor and its effective invariant, while actively computed (`src/mechanics/stress.jl`) for output, are not read back into the solve. `base_vertical`, `lateral`, and the two principal-stress eigenvalues are, in the current code, allocated but never written by anything in `src/`.

`effective` and `effective_depthaveraged` are two distinct quantities, not the same value at two resolutions: `effective` is DIVA's per-layer invariant ([robinson_comparison_2022](@citet) Eq. 13, including the vertical-shear terms) and feeds the column viscosity $\eta(z)$; `effective_depthaveraged` is the SSA invariant (Eq. 12, shear terms dropped) and feeds $\bar\eta$. Both are live.

$\bar u, \bar v$ are the prognostic unknowns SSA/DIVA solve for and directly build the flux $q = H\bar u$; $u_{\mathrm{b}}, v_{\mathrm{b}}$ feed the basal friction law — but only the vector components are currently written, not the `base_norm` magnitude. $u(z), v(z)$ are prognostic unknowns for Blatter–Pattyn; $w(z)$ is diagnosed from incompressibility (`verticalvelocity!`) but is itself a necessary input to the vertical-shear strain rate $\dot\varepsilon_{xz}, \dot\varepsilon_{yz}$ and hence to $\eta(z)$ — all three are live. All the velocity gradients are intermediates on the way to the strain-rate tensor and are live. Surface velocity ($u_{\mathrm{s}}, v_{\mathrm{s}}$, but again not `surface_norm`) is computed and kept up to date for observational comparison even though it is not read back into the solve; the separate 3D velocity-magnitude field `norm` is allocated but unused.

$F_1, F_2 = \int_b^s \mu^{-1}\left(\frac{s-z}{H}\right)^m dz$ ([robinson_comparison_2022](@citet) Eq. 15) are necessary for the DIVA momentum solve specifically, and are live (`viscosity_integrals!`) — they stay zero on any state that never runs DIVA. `rate_factor`/`rate_factor_depthaveraged` are, per the docstring at `src/api/state.jl:212`, prescribed inputs rather than derived from temperature — they are read (not computed) to drive the pseudo-transient solver's `GlenViscosityContinuation`, so "necessary" here means necessary for that specific solver path, not part of the temperature-to-viscosity chain (`ThermodynamicState` is currently an unconnected sibling; see below).

## [Thermodynamics](@id thermodynamics_state)

!!! note "Not yet wired to a solver"
    `AbstractEnergySolver` currently has no concrete subtype, and `TemperatureEnergyBalance`/`EnthalpyEnergyBalance` are marker structs with no methods — the heat equation is not yet implemented. Every field below is allocated in `ThermodynamicState` but not read or written anywhere else in `src/`, hence **Used** is ✗ throughout this table. **Necessary** instead reflects each field's intended role in the heat equation once it lands.

| Variable | Symbol | Dimension | Instances | Necessary | Used |
|---|---|---|---|---|---|
| ice temperature | $T$ | 2D/3D | `TemperatureState` | ✓ | ✗ |
| homologous temperature | $T'$ | 2D/3D | `TemperatureState` | ✓ | ✗ |
| surface temperature | $T_{\mathrm{s}}$ | 2D | `TemperatureState` | ✓ | ✗ |
| bedrock temperature | $T_{\mathrm{rock}}$ | 2D/3D | `TemperatureState` | ✓ | ✗ |
| pressure melting point | $T_{\mathrm{m}}(p)$ | 2D/3D | `TemperatureState` | ✓ | ✗ |
| ice enthalpy | $\mathcal{E}_{\mathrm{ice}}$ | 2D/3D | `EnthalpyState` | ✓ | ✗ |
| bedrock enthalpy | $\mathcal{E}_{\mathrm{rock}}$ | 2D/3D | `EnthalpyState` | ✓ | ✗ |
| ice water content | $\omega$ | 2D/3D | `ThermodynamicState` | ✓ | ✗ |
| internal strain heating | $\Phi$ | 2D/3D | `ThermodynamicState` | ✓ | ✗ |
| basal frictional heating | $\Phi_{\mathrm{b}}$ | 2D | `ThermodynamicState` | ✓ | ✗ |
| heat flux at ice base | $Q_{\mathrm{base}}$ | 2D | `ThermodynamicState` | ✓ | ✗ |
| heat flux into bedrock | $Q_{\mathrm{bedrock}}$ | 2D | `ThermodynamicState` | ✓ | ✗ |
| geothermal heat flux | $Q_{\mathrm{geo}}$ | 2D | `ThermodynamicState` | ✓ | ✗ |
| specific heat capacity of ice | $c_p$ | 2D/3D | `ThermodynamicState` | ✓ | ✗ |
| thermal conductivity of ice | $k$ | 2D/3D | `ThermodynamicState` | ✓ | ✗ |
| strain-heating tendency | $\partial \Phi/\partial t$ | 2D/3D | `ThermodynamicState` | ✗ | ✗ |
| basal water-layer thickness | $H_{\mathrm{w}}$ | 2D | `ThermodynamicState` | ✗ | ✗ |
| water-layer thickness tendency | $\partial H_{\mathrm{w}}/\partial t$ | 2D | `ThermodynamicState` | ✗ | ✗ |
| cold-temperate transition surface | $z_{\mathrm{CTS}}$ | 2D | `ThermodynamicState` | ✗ | ✗ |

$T$ is one of the model's four main prognostic variables. $T'$ (temperature relative to the pressure melting point) is what the rate factor $A(T')$ actually consumes ([Material](@ref material)), so it sits directly on the temperature-to-viscosity path. $T_{\mathrm{s}}$ and $T_{\mathrm{rock}}$ are the surface and basal boundary conditions the heat equation needs; $T_{\mathrm{m}}(p)$ caps $T$ and defines $T'$. The enthalpy formulation is an alternative prognostic representation of temperature (under `EnthalpyEnergyBalance`), playing the same role as `ice`/`rock` in `TemperatureState`.

$\omega$, $\Phi$, $\Phi_{\mathrm{b}}$, the three heat fluxes, and $c_p$/$k$ are source terms and coefficients the heat equation needs. $H_{\mathrm{w}}$ would become necessary if/when it feeds a hydrology-dependent effective-pressure or friction law; today it (and its tendency) are bookkeeping fields. $z_{\mathrm{CTS}}$ marks the polythermal cold/temperate interface for output and is not itself consumed elsewhere.

## Material

| Variable | Symbol | Dimension | Instances | Necessary | Used |
|---|---|---|---|---|---|
| column viscosity | $\eta(z)$ | 2D/3D | `MaterialState`, `MechanicMaterialState` | ✓ | ✓ |
| depth-averaged viscosity | $\bar\eta$ | 2D | `MaterialState`, `MechanicMaterialState` | ✓ | ✓ |
| depth-integrated viscosity | $\int\eta\,dz$ | 2D | `MaterialState` | ✗ | ✗ |

Viscosity is the clearest example of a variable duplicated by design: `MaterialState.eta_ice` is the material component's own copy (its own grid), currently feeding the diagnostic stress-tensor calculation in `src/mechanics/stress.jl`, while `MechanicMaterialState.viscosity`/`viscosity_depthaveraged` is the copy the momentum solver reads and updates directly — the two are not automatically kept in sync. Of the two `MaterialState` copies, only the column one (`eta_ice`) is actually used; `MaterialState.eta_depth_averaged` is currently unused (only `MechanicMaterialState.viscosity_depthaveraged` is live), which is why the depth-averaged row is marked used overall — at least one of its instances is. `eta_depth_integrated` ($\int_b^s\eta\,dz$, distinct from the depth-*average* $\bar\eta = H^{-1}\int_b^s\eta\,dz$) is allocated but not read anywhere in `src/`.

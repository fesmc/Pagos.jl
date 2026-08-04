"""
Seconds per Julian year (365.25 d), the conversion between the published SI convention
and Pagos's internal one. See [`Constants`](@ref) for what "internal" means here.
"""
const SECONDS_PER_YEAR = 31_557_600.0

"""
$(TYPEDSIGNATURES)

Physical constants that do not depend on any modelling choice.

!!! note "Base units are metre, **year**, pascal"
    Pagos computes in `(m, yr, Pa)`, not full SI: velocities are `m yr⁻¹`, strain rates
    `yr⁻¹`, viscosities `Pa yr`, rate factors `Pa⁻ⁿ yr⁻¹`, and timesteps `yr`.

    Taking `Pa` as a base unit rather than deriving it from `kg` is what keeps this
    cheap: stress is then unaffected by the time unit, so `density_*` and `gravity`
    below keep their ordinary SI values. They are only ever used as the product `ρg`
    (a specific weight, `Pa m⁻¹`, which carries no time) or as dimensionless density
    ratios in the flotation criterion — never separately in a way that would expose
    `kg m⁻¹ yr⁻²`. The same applies to every other field here: activation energies,
    heat capacities and conductivities are per-kelvin/per-kilogram/per-mole quantities
    with no time dimension.

    What *does* carry the time unit lives in the rheology — see
    [`AbstractRateFactor`](@ref) and the `time_unit` constructor argument for passing
    values published in `s⁻¹`.

Override selectively via keyword arguments:

```julia
cst = Constants{Float64}()                          # standard defaults
cst = Constants{Float32}(density_ice = Float32(910))
```
"""
@kwdef struct Constants{T<:AbstractFloat}
    # Densities (kg m⁻³)
    density_ice::T        = 910.0
    density_seawater::T   = 1028.0
    density_freshwater::T = 1000.0
    density_rock::T       = 3300.0

    # Mechanics
    gravity::T            = 9.81         # gravitational acceleration (m s⁻²)

    # Thermodynamics
    melting_temperature::T              = 273.15   # reference melting point at standard pressure (K)
    gas_constant::T                     = 8.314    # universal gas constant (J mol⁻¹ K⁻¹)
    latent_heat_fusion_ice::T           = 3.34e5   # latent heat of fusion of ice (J kg⁻¹)
    specific_heat_capacity_ice::T       = 2009.0   # (J kg⁻¹ K⁻¹)
    specific_heat_capacity_seawater::T  = 3974.0   # (J kg⁻¹ K⁻¹)
    specific_heat_capacity_rock::T      = 1000.0   # (J kg⁻¹ K⁻¹)
    thermal_conductivity_ice::T         = 2.1      # (W m⁻¹ K⁻¹)
    thermal_conductivity_rock::T        = 3.0      # (W m⁻¹ K⁻¹)
    clausius_clapeyron_slope::T         = 9.8e-8   # pressure-melting slope (K Pa⁻¹)

    # Unit conversion
    seconds_per_year::T                 = SECONDS_PER_YEAR   # (s yr⁻¹)
end

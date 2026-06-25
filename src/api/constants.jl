"""
$(TYPEDSIGNATURES)

Physical constants that do not depend on any modelling choice.
All values are in SI units. Override selectively via keyword arguments:

```julia
cst = Constants{Float64}()                          # standard defaults
cst = Constants{Float32}(density_ice = Float32(910))
```
"""
@kwdef struct Constants{T<:AbstractFloat}
    # Densities (kg m⁻³)
    density_ice::T        = T(910.0)
    density_seawater::T   = T(1028.0)
    density_freshwater::T = T(1000.0)
    density_rock::T       = T(3300.0)

    # Mechanics
    gravity::T            = T(9.81)         # gravitational acceleration (m s⁻²)

    # Thermodynamics
    melting_temperature::T              = T(273.15)   # reference melting point at standard pressure (K)
    gas_constant::T                     = T(8.314)    # universal gas constant (J mol⁻¹ K⁻¹)
    latent_heat_fusion_ice::T           = T(3.34e5)   # latent heat of fusion of ice (J kg⁻¹)
    specific_heat_capacity_ice::T       = T(2009.0)   # (J kg⁻¹ K⁻¹)
    specific_heat_capacity_seawater::T  = T(3974.0)   # (J kg⁻¹ K⁻¹)
    specific_heat_capacity_rock::T      = T(1000.0)   # (J kg⁻¹ K⁻¹)
    thermal_conductivity_ice::T         = T(2.1)      # (W m⁻¹ K⁻¹)
    thermal_conductivity_rock::T        = T(3.0)      # (W m⁻¹ K⁻¹)
    clausius_clapeyron_slope::T         = T(9.8e-8)   # pressure-melting slope (K Pa⁻¹)

    # Unit conversion
    seconds_per_year::T                 = T(31_557_600.0)   # (s yr⁻¹)
end

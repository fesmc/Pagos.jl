"""
$(TYPEDSIGNATURES)

Physical constants used in Pagos, with default values.

# Fields
 - `ρ_ice`: Density of ice (kg m⁻³)
 - `ρ_water`: Density of ocean water (kg m⁻³)
 - `ρ_freshwater`: Density of freshwater (kg m⁻³)
 - `g`: Gravitational acceleration (m s⁻²)
 - `T0`: Reference melting temperature (K)
 - `R`: Universal gas constant (J mol⁻¹ K⁻¹)
 - `L_ice`: Latent heat of fusion of ice (J kg⁻¹)
 - `c_ice`: Specific heat capacity of ice (J kg⁻¹ K⁻¹)
 - `k_ice`: Thermal conductivity of ice (W m⁻¹ K⁻¹)
 - `β_cc`: Clausius-Clapeyron slope (K Pa⁻¹)
 - `c_ocean`: Specific heat capacity of seawater (J kg⁻¹ K⁻¹)
 - `spy`: Seconds per year (s yr⁻¹)
"""
@kwdef struct PhysicalConstants{T}
    # Densities
    ρ_ice::T        = 910.0
    ρ_water::T      = 1028.0
    ρ_freshwater::T = 1000.0

    # Mechanics
    g::T            = 9.81

    # Thermodynamics
    T0::T           = 273.15
    R::T            = 8.314
    L_ice::T        = 3.34e5
    c_ice::T        = 2009.0
    k_ice::T        = 2.1
    β_cc::T         = 9.8e-8
    c_ocean::T      = 3974.0

    # Unit conversion
    spy::T          = 31_557_600.0
end

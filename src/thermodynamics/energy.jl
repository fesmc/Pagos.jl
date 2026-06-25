abstract type AbstractEnergyBalance end

"""
$(TYPEDSIGNATURES)

Enforce the ice thermodynamics to be constant (no energy balance).
"""
struct NoEnergyBalance <: AbstractEnergyBalance end

"""
$(TYPEDSIGNATURES)

Model the ice thermodynamics via temperature-based energy balance.
"""
struct TemperatureEnergyBalance <: AbstractEnergyBalance end

"""
$(TYPEDSIGNATURES)

Model the ice thermodynamics via enthalpy-based energy balance.
"""
struct EnthalpyEnergyBalance <: AbstractEnergyBalance end
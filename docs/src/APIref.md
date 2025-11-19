# API reference

## Dynamics

### Effective Pressure

```@docs
AbstractEffectivePressure
ConstantEffectivePressure
OverburdenEffectivePressure
effective_pressure
```

### Friction

```@docs
AbstractBasalFriction
ConstantBetaBasalFriction
LinearBetaBasalFriction
PseudoPlasticPowerBasalFriction
RegularizedCoulombBasalFriction
basal_shear_stress
basal_shear_stress!
```

## Calving

```@docs
AbstractCalving
LipscombCalving
LevermannCalving
calving_rate
```

## Material

### Pressure melting point

```@docs
AbstractPressureMeltingPoint
LinearPressureMeltingPoint
```

### Rate factor

```@docs
AbstractRateFactor
ConstantRateFactor
ArrheniusRateFactor
SmithMorlandRateFactor
rate_factor
rate_factor!
```

### Creep function

```@docs
AbstractCreepFunction
GlenNyeCreepFunction
RegularizedGlenNyeCreepFunction
SmithMorlandCreepFunction
creep_function
creep_function!
```

### Flow law

```@docs
AbstractFlowLaw
ConstantViscosityFlowLaw
RateCreepFlowLaw
GlenNyeFlowLaw
RegularizedGlenNyeFlowLaw
SmithMorlandFlowLaw
```

## Plots

```@docs
plot_rate_factor
plot_melting_point
plot_ice_viscosity
plot_basal_shear_stress
```
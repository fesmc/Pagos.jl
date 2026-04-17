# API reference

## Dynamics

### Ice Dynamics

```@docs
AbstractDynamics
SIADynamics
SSADynamics
HybridDynamics
L1L2Dynamics
DIVADynamics
BlatterPattynDynamics
StokesDynamics
velocity
velocity!
```

### Solvers

```@docs
AbstractDynamicsSolver
MatrixDynamicsSolver2D
MatrixDynamicsSolver3D
PseudoTransientDynamicsSolver2D
PseudoTransientDynamicsSolver3D
NeuralDynamicsSolver2D
NeuralDynamicsSolver3D
```

### Effective Pressure

```@docs
AbstractEffectivePressure
ConstantEffectivePressure
OverburdenEffectivePressure
effective_pressure
effective_pressure!
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

## Topography

### Sigma transform

```@docs
AbstractSigmaTransform
LinearSigmaTransform
ExponentialSigmaTransform
sigma_transform
```

### Calving

```@docs
AbstractCalving
LipscombCalving
LevermannCalving
calving_rate
calving_rate!
```

## Material

### Pressure melting point

```@docs
AbstractPressureMeltingPoint
LinearPressureMeltingPoint
pressure_melting_point
pressure_melting_point!
relative_temperature
relative_temperature!
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

## Thermodynamics

## Boundary conditions


## Plots

```@docs
plot_rate_factor
plot_melting_point
plot_ice_viscosity
plot_basal_shear_stress
```
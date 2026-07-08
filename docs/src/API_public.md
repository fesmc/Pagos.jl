# Public API

## Dynamics

### Ice Dynamics

```@docs
AbstractMomentumBalance
SIAMomentumBalance
SSAMomentumBalance
HybridMomentumBalance
L1L2MomentumBalance
DIVAMomentumBalance
BlatterPattynMomentumBalance
StokesMomentumBalance
velocity
velocity!
```

### Solvers

```@docs
AbstractMomentumSolver
MatrixMomentumSolver2D
MatrixMomentumSolver3D
PseudoTransientMomentumSolver2D
PseudoTransientMomentumSolver3D
NeuralMomentumSolver2D
NeuralMomentumSolver3D
```

### Effective Pressure

```@docs
AbstractEffectivePressure
PrescribedEffectivePressure
OverburdenEffectivePressure
LeguyEffectivePressure
TillEffectivePressure
effective_pressure
effective_pressure!
```

### Friction

```@docs
AbstractBasalBeta
PrescribedBasalBeta
PseudoPlasticPowerBasalBeta
CoulombBasalBeta
basal_shear_stress
basal_shear_stress!
```

## Topography

### Sigma transform

```@docs
AbstractSigmaTransform
PowerSigmaTransform
LinearSigmaTransform
QuadraticSigmaTransform
VerticalLayering
CorrectedVerticalLayering
sigma
get_ζ_aa
get_ζ_ac
```

### Calving

```@docs
AbstractCalving
PrescribedCalving
RelaxedCalving
ThicknessCalving
FlotationCalving
LipscombCalving
LevermannCalving
CrawfordCalving
BassisCalving
EigenCalving
PollardDeContoCalving
calving_rate
calving_rate!
```

## Material

### Pressure melting point

```@docs
AbstractPressureMeltingPoint
PrescribedPressureMeltingPoint
LinearPressureMeltingPoint
LinearSalinityPressureMeltingPoint
pressure_melting_point
pressure_melting_point!
relative_temperature
relative_temperature!
thermal_forcing
thermal_forcing!
```

### Rate factor

```@docs
AbstractRateFactor
PrescribedRateFactor
ArrheniusRateFactor
SmithMorlandRateFactor
HookeRateFactor
LliboutryDuvalRateFactor
FanLowStrainGSIRateFactor
FanLowStrainGSS1RateFactor
FanLowStrainGSS2RateFactor
FanHighStrainGSIRateFactor
rate_factor
rate_factor!
```

### Creep function

```@docs
AbstractCreep
GlenNyeCreep
SmithMorlandCreep
PettitWaddingtonCreep
GoldsbyKohlstedtCreep
FanLowStrainCreep
creep
creep!
```

### Flow law

```@docs
AbstractFlowLaw
PrescribedViscosityFlowLaw
RateCreepFlowLaw
GlenNyeFlowLaw
SmithMorlandFlowLaw
FanLowStrainFlowLaw
FanHighStrainFlowLaw
viscosity
viscosity!
```

### Anisotropy

```@docs
AbstractAnisotropy
EnhancementAnisotropy
CAFFEAnisotropy
anisotropy!
enhancement_factor
deformability
square_tangential_invariant
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

## Numerics

### Derivatives

```@docs
∂x!
∂y!
∂x₃!
∂x₁₂!
∂x₁
∂x₂
∂x₃
∂x₁₂
```

## Utilities

### Indexing

```@docs
AbstractIndexing
StrictIndexing
FlatIndexing
ReflectiveIndexing
PeriodicIndexing
index
stencil_fd
stencil
```

### Active cells map

```@docs
ActiveCellsMap
apply!
```

### Staggered grids

```@docs
StaggeredGrids
StaggeredGrids.RectilinearGrid
StaggeredGrids.Field
StaggeredGrids.interpolate
StaggeredGrids.compute
StaggeredGrids.compute!
StaggeredGrids.@at
```
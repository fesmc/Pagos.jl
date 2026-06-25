module Pagos

using DocStringExtensions
using Downloads
using KernelAbstractions
using LinearAlgebra
using LinearSolve           # TODO: externalize
using LoopVectorization
using Random
using StatsBase
using SparseArrays
using Tullio

###########################################################
# Structs
###########################################################

include("api/domain.jl")
include("api/state.jl")
include("api/constants.jl")
include("api/icesheet.jl")
include("api/simulation.jl")
include("api/step.jl")
# TODO: remove Params — physics constants moved to Constants, solver params live in solver structs
# include("api/params.jl")
# TODO: remove Options — all fields already live in PseudoTransientSolver
# include("api/options.jl")
export AbstractGrid, RegularGrid, CommonGrid
export TopographyState, DynamicsState, ThermodynamicsState, MaterialState
export Constants, IceSheet, Simulation

###########################################################
# Utils
###########################################################

include("utils/indices.jl")
export AbstractIndexing, StrictIndexing, FlatIndexing
export ReflectiveIndexing, PeriodicIndexing
export index, stencil_fd, stencil

include("utils/math.jl")
include("utils/debug.jl")
include("utils/mask.jl")
export ActiveCellsMap, active_indices!, apply!

include("utils/integrators.jl")
export AbstractIntegrationMethod, Euler, RungeKutta4, BogackiShampine32, Tsitouras54, RKL2
export SSPRK33, SSPRK43, ssp_coefficient
export Integrator, step!, substep!, estimate_spectral_radius!

###########################################################
# Topography
###########################################################

include("topography/sigmatransform.jl")
export AbstractSigmaTransform, PowerSigmaTransform
export LinearSigmaTransform, QuadraticSigmaTransform
export VerticalLayering, CorrectedVerticalLayering
export ζ_aa, ζ_ac, sigma
export get_ζ_aa, get_ζ_ac

include("topography/calving.jl")
export AbstractCalving, AbstractGroundedCalving, AbstractFloatingCalving
export PrescribedCalving, RelaxedCalving, ThicknessCalving
export LipscombCalving, LevermannCalving, CrawfordCalving, BassisCalving
export EigenCalving, FlotationCalving, PollardDeContoCalving
export calving_rate, calving_rate!

# TODO: update to new API (State → TopographyState/DynamicsState, Domain → RegularGrid)
# include("topography/advection.jl")

###########################################################
# Dynamics
###########################################################

include("mechanics/velocities.jl")
export AbstractDynamics
export SIADynamics, SSADynamics, SIASSADynamics
export DIVADynamics, BlatterPattynDynamics, StokesDynamics
export InertialDIVADynamics, InertialSIASSADynamics
export populate_vectors!, velocity!
export AbstractDynamicsSolver, LinearDynamicsSolver2D

include("legacy/performance/velocities.jl")
export LegacyLinearDynamicsSolver2D

include("mechanics/effective_pressure.jl")
export AbstractEffectivePressure, PrescribedEffectivePressure
export OverburdenEffectivePressure
export LeguyEffectivePressure, TillEffectivePressure
export effective_pressure, effective_pressure!

include("mechanics/basal_friction.jl")
export AbstractFriction, BasalFriction
export AbstractBasalBeta, PrescribedBasalBeta
export PseudoPlasticPowerBasalBeta, CoulombBasalBeta
export basal_shear_stress, basal_shear_stress!

include("mechanics/velocities3D.jl")
export aggregate_viscosity_integral!, aggregated_viscosity_integral!
export velocities3D!, surface_velocity!, depthaveraged_velocity!
export basal_velocity_from_surface_velocity!, depthavg_velocity!

# TODO: update to new API (State → DynamicsState, Domain → RegularGrid)
# include("mechanics/pseudotransient.jl")
# export pseudo_transient!, pseudo_dotvel!, pseudo_dt, pseudo_vel!, dotvel!

include("mechanics/inertial.jl")
export inertial_velocity!

# include("mechanics/friction/plastic.jl")
# include("mechanics/friction/stagger.jl")

###########################################################
# Material
###########################################################

include("material/anisotropy.jl")
export AbstractAnisotropy, EnhancementAnisotropy, CAFFEAnisotropy
export anisotropy!, enhancement_factor
export deformability, square_tangential_invariant

include("material/rate_factor.jl")
export AbstractRateFactor, PrescribedRateFactor, ArrheniusRateFactor
export SmithMorlandRateFactor, HookeRateFactor, LliboutryDuvalRateFactor
export FanLowStrainGSIRateFactor, FanLowStrainGSS1RateFactor, FanLowStrainGSS2RateFactor
export FanHighStrainGSIRateFactor
export rate_factor, rate_factor!

include("material/creep.jl")
export AbstractCreep, GlenNyeCreep, SmithMorlandCreep
export PettitWaddingtonCreep, GoldsbyKohlstedtCreep
export FanLowStrainCreep
export creep, creep!

include("material/flow_law.jl")
export AbstractFlowLaw, PrescribedViscosityFlowLaw, RateCreepFlowLaw
export GlenNyeFlowLaw, SmithMorlandFlowLaw
export FanLowStrainFlowLaw, FanHighStrainFlowLaw
export viscosity, viscosity!

include("material/pressure_melting_point.jl")
export AbstractPressureMeltingPoint, LinearPressureMeltingPoint
export PrescribedPressureMeltingPoint, LinearSalinityPressureMeltingPoint
export pressure_melting_point, pressure_melting_point!
export relative_temperature, relative_temperature!
export thermal_forcing, thermal_forcing!

include("material/stress.jl")
export shearstress!, basalstress!, drivingstress!

include("material/strainrate.jl")
export scaledstrainrate!, velocitygradients!

###########################################################
# Numerics
###########################################################

include("numerics/differences.jl")
include("numerics/picard.jl")
# include("numerics/pseudotransient.jl")
# export stagger_beta!
# export pseudo_dotvel!, pseudo_vel!, pseudo_transient!
export ∂x₁, ∂x₂, ∂x₃, ∂x₁!, ∂x₂!, ∂x₃!
export ∂x₁₂, ∂x₁₂!

###########################################################
# Extensions
###########################################################

include("ext/plots.jl")

end
module Pagos

using DocStringExtensions
using Downloads
using KernelAbstractions
using LinearAlgebra
using LinearSolve           # TODO: externalize
using StatsBase
using SparseArrays

###########################################################
# Structs
###########################################################

include("structs/domain.jl")
include("structs/state.jl")
include("structs/params.jl")
include("structs/options.jl")
# include("structs/tools.jl")
include("structs/icesheet.jl")
export Domain, State, Params, Options, IceSheet

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

include("topography/advection.jl")

###########################################################
# Dynamics
###########################################################

include("dynamics/effective_pressure.jl")
export AbstractEffectivePressure, PrescribedEffectivePressure
export OverburdenEffectivePressure
export LeguyEffectivePressure, TillEffectivePressure
export effective_pressure, effective_pressure!

include("dynamics/basal_friction.jl")
export AbstractBasalBeta, PrescribedBasalBeta
export PseudoPlasticPowerBasalBeta, CoulombBasalBeta
export basal_shear_stress, basal_shear_stress!

include("dynamics/velocities3D.jl")
export aggregate_viscosity_integral!, aggregated_viscosity_integral!
export velocities3D!, surface_velocity!, depthaveraged_velocity!
export basal_velocity_from_surface_velocity!, depthavg_velocity!

include("dynamics/pseudotransient.jl")
export pseudo_transient!

include("dynamics/inertial.jl")
export inertial_velocity!

include("dynamics/friction/plastic.jl")
include("dynamics/friction/stagger.jl")

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
include("material/strainrate.jl")

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
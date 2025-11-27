module Pagos

using DocStringExtensions
using Downloads
using KernelAbstractions
using LinearAlgebra
using LinearSolve           # TODO: externalize
using StatsBase
using SparseArrays
# using CUDA

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
include("utils/math.jl")
include("utils/debug.jl")

# include("helpers/indices.jl")
# include("helpers/staggering.jl")

###########################################################
# Topography
###########################################################

include("topography/sigmatransform.jl")
export AbstractSigmaTransform, PowerSigmaTransform, ArctanSigmaTransform
export VerticalLayering, CorrectedVerticalLayering
export ζ_aa, ζ_ac, sigma

include("topography/calving.jl")
export AbstractCalving, RelaxedCalving, ThicknessCalving, LipscombCalving
export LevermannCalving, CrawfordCalving, BassisCalving, BedStddevCalving
export calving_rate, calving_rate!

###########################################################
# Dynamics
###########################################################

include("dynamics/effective_pressure.jl")
export AbstractEffectivePressure, ConstantEffectivePressure
export OverburdenEffectivePressure
export LeguyEffectivePressure, TillEffectivePressure
export effective_pressure, effective_pressure!

include("dynamics/basal_friction.jl")
export AbstractBasalFriction, ConstantBetaBasalFriction
export LinearBetaBasalFriction
export PseudoPlasticPowerBasalFriction, RegularizedCoulombBasalFriction
export basal_shear_stress, basal_shear_stress!

include("dynamics/velocities3D.jl")
export aggregate_viscosity_integral!, aggregated_viscosity_integral!
export velocities3D!, surface_velocity!, depthaveraged_velocity!
export basal_velocity_from_surface_velocity!, depthavg_velocity!

include("dynamics/friction/plastic.jl")
include("dynamics/friction/stagger.jl")
include("dynamics/advection.jl")

###########################################################
# Material
###########################################################

include("material/rate_factor.jl")
export AbstractRateFactor, ConstantRateFactor, ArrheniusRateFactor
export SmithMorlandRateFactor
export rate_factor, rate_factor!

include("material/creep_function.jl")
export AbstractCreepFunction, GlenNyeCreepFunction, RegularizedGlenNyeCreepFunction
export SmithMorlandCreepFunction
export creep_function, creep_function!

include("material/flow_law.jl")
export AbstractFlowLaw, ConstantViscosityFlowLaw, RateCreepFlowLaw
export GlenNyeFlowLaw, RegularizedGlenNyeFlowLaw, SmithMorlandFlowLaw
export viscosity, viscosity!

include("material/pressure_melting_point.jl")
export AbstractPressureMeltingPoint, LinearPressureMeltingPoint
export melting_point, melting_point!
export relative_temperature, relative_temperature!

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
# export delx, dely, advect!

###########################################################
# Extensions
###########################################################

include("ext/plots.jl")

end
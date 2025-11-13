module Pagos

using DocStringExtensions
using Downloads
using LinearAlgebra
using StatsBase
# using LinearSolve
# using SparseArrays
# using CUDA

include("structs/domain.jl")
include("structs/state.jl")
include("structs/params.jl")
include("structs/options.jl")
# include("structs/tools.jl")
include("structs/icesheet.jl")
export Domain, State, Params, Options, IceSheet

include("helpers/indices.jl")
include("helpers/staggering.jl")
include("helpers/debug.jl")
include("helpers/sigmatransform.jl")
export exponential_vertical_layers, sigma

###########################################################
# Dynamics
###########################################################

include("dynamics/effective_pressure.jl")
export AbstractEffectivePressure, ConstantEffectivePressure
export OverburdenEffectivePressure
export LeguyEffectivePressure, TillEffectivePressure
export effective_pressure

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
export rate_factor!, rate_factor, matrix_rate_factor

include("material/creep_function.jl")
export AbstractCreepFunction, GlenNyeCreepFunction, RegularizedGlenNyeCreepFunction
export SmithMorlandCreepFunction
export creep_function!, creep_function

include("material/flow_law.jl")
export AbstractFlowLaw, ConstantViscosityFlowLaw, RateCreepFlowLaw
export GlenNyeFlowLaw, RegularizedGlenNyeFlowLaw, SmithMorlandFlowLaw
export viscosity



include("material/pressure_melting_point.jl")
export AbstractPressureMeltingPoint, LinearPressureMeltingPoint
export melting_point, melting_point!
export get_relative_temperature, update_relative_temperature!

include("material/stress.jl")
include("material/strainrate.jl")

###########################################################
# Numerics
###########################################################

include("numerics/differences.jl")
include("numerics/picard.jl")
include("numerics/pseudotransient.jl")
export stagger_beta!
export pseudo_dotvel!, pseudo_vel!, pseudo_transient!
export delx, dely, advect!

###########################################################
# Extensions
###########################################################

include("ext/plots.jl")

end
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
include("utils/math.jl")
include("utils/debug.jl")
include("utils/mask.jl")
export ActiveCellsMap, active_indices!, apply!

include("helpers/indices.jl")
# include("helpers/staggering.jl")

###########################################################
# Topography
###########################################################

include("topography/sigmatransform.jl")
export AbstractSigmaTransform, PowerSigmaTransform
export ArctanSigmaTransform, LinearSigmaTransform
export QuadraticSigmaTransform
export VerticalLayering, CorrectedVerticalLayering
export ζ_aa, ζ_ac, sigma

include("topography/calving.jl")
export AbstractCalving, AbstractGroundedCalving, AbstractFloatingCalving
export ConstantCalving, RelaxedCalving, ThicknessCalving
export LipscombCalving, LevermannCalving, CrawfordCalving, BassisCalving
export EigenCalving, FlotationCalving, PollardDeContoCalving
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
export AbstractBasalBeta, ConstantBasalBeta
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
include("dynamics/advection.jl")

###########################################################
# Material
###########################################################

include("material/rate_factor.jl")
export AbstractRateFactor, ConstantRateFactor, ArrheniusRateFactor
export SmithMorlandRateFactor, HookeRateFactor, LliboutryDuvalRateFactor
export rate_factor, rate_factor!

include("material/creep.jl")
export AbstractCreep, GlenNyeCreep, RegularizedGlenNyeCreep
export SmithMorlandCreep, PettitWaddingtonCreep, GoldsbyKohlstedtCreep
export creep, creep!

include("material/flow_law.jl")
export AbstractFlowLaw, ConstantViscosityFlowLaw, RateCreepFlowLaw
export GlenNyeFlowLaw, RegularizedGlenNyeFlowLaw, SmithMorlandFlowLaw
export viscosity, viscosity!

include("material/pressure_melting_point.jl")
export AbstractPressureMeltingPoint, LinearPressureMeltingPoint
export ConstantPressureMeltingPoint, LinearSalinityPressureMeltingPoint
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
export delx1, delx2, delx1!, delx2!

###########################################################
# Extensions
###########################################################

include("ext/plots.jl")

end
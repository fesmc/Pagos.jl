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

include("dynamics/friction/plastic.jl")
include("dynamics/friction/stagger.jl")
include("dynamics/advection.jl")

include("dynamics/velocities3D.jl")
export aggregate_viscosity_integral!, aggregated_viscosity_integral!
export velocities3D!, surface_velocity!, depthaveraged_velocity!
export basal_velocity_from_surface_velocity!, depthavg_velocity!

include("material/stress.jl")
include("material/strainrate.jl")

include("numerics/differences.jl")
include("numerics/picard.jl")
include("numerics/pseudotransient.jl")
export stagger_beta!
export pseudo_dotvel!, pseudo_vel!, pseudo_transient!

export delx, dely, advect!

include("material/viscosity.jl")
export AbstractViscosity, GlenNyeViscosity, RegularizedGlenNyeViscosity, ConstantViscosity
export get_viscosity

include("material/rate_factor.jl")
export AbstractRateFactor, ConstantRateFactor, ArrheniusRateFactor, SmithMorlandRateFactor
export update_rate_factor!, get_rate_factor

include("material/pressure_melting_point.jl")
export AbstractPressureMeltingPoint, LinearPressureMeltingPoint
export get_melting_point, update_melting_point!
export get_relative_temperature, update_relative_temperature!


include("plots.jl")

end
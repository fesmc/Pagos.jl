module Pagos

using LinearAlgebra
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

include("thermodynamics/viscosity.jl")
export delx, dely, advect!

end
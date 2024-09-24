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

include("helpers/indices.jl")
include("helpers/staggering.jl")
include("helpers/debug.jl")

include("dynamics/velocities3D.jl")
include("dynamics/friction/plastic.jl")
include("dynamics/friction/stagger.jl")
include("dynamics/advection.jl")

include("material/stress.jl")
include("material/strainrate.jl")

include("numerics/differences.jl")
include("numerics/picard.jl")
include("numerics/pseudotransient.jl")

include("thermodynamics/viscosity.jl")

export Domain, State, Params, Options, IceSheet
export stagger_beta!
export pseudo_dotvel!, pseudo_vel!, pseudo_transient!
export delx, dely, advect!
export aggregate_viscosity_integral!, velocities3D!
export surface_velocity!, depthaveraged_velocity!

end
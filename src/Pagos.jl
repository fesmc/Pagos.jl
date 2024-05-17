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

include("mechanics/velocity.jl")
include("mechanics/stress.jl")
include("mechanics/strainrate.jl")
include("mechanics/friction/plastic.jl")
include("mechanics/friction/stagger.jl")

include("numerics/differences.jl")
include("numerics/picard.jl")
include("numerics/pseudotransient.jl")

include("thermodynamics/viscosity.jl")

export Domain, State, Params, Options, IceSheet

export pseudo_dotvel!, pseudo_vel!, pseudo_transient!

end
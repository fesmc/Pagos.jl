abstract type AbstractGrainSize end

struct PrescribedGrainSize{M} <: AbstractGrainSize
    grain_size::M   # <: AbstractArray or Real
end
abstract type AbstractConvolution end

struct PrecomputedFFTConvolution{M} <: AbstractConvolution
    kernel::M   # <: AbstractArray
end
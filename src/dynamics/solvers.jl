"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the dynamics solver via [`veolicity`](@ref).

# Available subtypes:
"""
abstract type AbstractDynamicsSolver end

struct MatrixDynamicsSolver end
struct OptimDynamicsSolver end
struct WaveletDynamicsSolver end
struct PseudoTransientDynamicsSolver end
struct ConvolutionalDynamicsSolver end
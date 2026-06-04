"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the dynamics solver via [`velocity`](@ref).

# Available subtypes:
"""
abstract type AbstractDynamicsSolver end

"""
$(TYPEDSIGNATURES)

TODO: Solve the ice dynamics via energy minimization.
"""
struct OptimDynamicsSolver <: AbstractDynamicsSolver
end

"""
$(TYPEDSIGNATURES)

TODO: Solve the ice dynamics via an iterative linear solver (e.g., CG, GMRES).
"""
struct IterativeDynamicsSolver <: AbstractDynamicsSolver
end

"""
$(TYPEDSIGNATURES)

TODO: Solve the ice dynamics via wavelet methods.
"""
struct WaveletDynamicsSolver <: AbstractDynamicsSolver
end

"""
$(TYPEDSIGNATURES)

TODO: Solve the ice dynamics via convolutional neural network (IGM style).
"""
struct ConvolutionalDynamicsSolver <: AbstractDynamicsSolver
end

"""
$(TYPEDSIGNATURES)

TODO: Solve the ice dynamics via a transient solver (e.g., explicit time-stepping).
"""
struct TransientDynamicsSolver <: AbstractDynamicsSolver
end

"""
$(TYPEDSIGNATURES)

TODO: Solve the ice dynamics via a pseudo-transient solver.
"""
struct PseudoTransientSolver <: AbstractDynamicsSolver
end

"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via a direct linear solver (e.g., sparse LU factorization).

# Improvements over `LegacyLinearDynamicsSolver2D`:
 1. dynamics is a type parameter → dispatch on loop1!/loop2! is fully static, no need to thread a runtime `dynamics` argument through populate_vectors! layers.
 2. SparseMatrixCSC pre-allocated at construction (fixed sparsity pattern). The hot path writes directly to A.nzval via a precomputed COO→nzval index map, eliminating the sparse(Ai, Aj, Av) allocation on every solve.
 3. Single AI type parameter (i_idx and j_idx are always the same kind).
 4. VT/MT/PI type parameters for vector/matrix/perm arrays so the struct can hold GPU arrays (CuVector, CuSparseMatrix) without code changes. The populate_vectors! kernels are written with KernelAbstractions and run on whichever backend owns lsd.u.
"""
struct LinearDynamicsSolver2D{
    DYN <: AbstractDynamics,
    RP  <: ResolutionParameters,
    VT  <: AbstractVector,        # float vector type (u, u0, b)
    MT,                            # sparse matrix type (SparseMatrixCSC or CuSparseMatrix)
    PI  <: AbstractVector{Int},   # perm index vector type
    AI,
} <: AbstractDynamicsSolver
    dynamics::DYN
    resolution_params::RP
    u::VT
    u0::VT
    b::VT
    A::MT
    perm::PI                      # COO fill order → A.nzval index
    i_idx::AI
    j_idx::AI
    solver_cache::Ref{Any}        # holds a cached direct solver (Nothing or backend-specific)
end
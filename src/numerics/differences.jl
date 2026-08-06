"""
$(TYPEDSIGNATURES)

In-place partial derivative of `u` along the first dimension. Uses central differences in
the interior; boundary behaviour is controlled by `idx` (default: `FlatIndexing`,
one-sided differences). GPU-compatible and asynchronous (see
[Asynchronous kernel launches](@ref)).
"""
function ∂x!(du, u, dx, idx::AbstractIndexing = FlatIndexing(1, size(u, 1)))
    backend = get_backend(u)
    if backend isa KernelAbstractions.CPU
        # CPU: 1D column kernel — stencil_fd(i, idx) is loop-invariant in the inner
        # i-loop, letting LLVM use a scalar reciprocal and SIMD multiplications
        # instead of per-lane vector divisions.  GPU needs the 2D kernel for
        # full thread occupancy, so the fast path is CPU-only.
        kernel! = _∂x₁_col!(backend)
        kernel!(du, u, dx, idx; ndrange = size(u, 2))
    else
        kernel! = _∂x!(backend)
        kernel!(du, u, dx, idx; ndrange = size(u))
    end
    return nothing
end

# Backward-compatible dispatch: called with an integer grid size instead of an indexing.
∂x!(du, u, dx, ::Integer) = ∂x!(du, u, dx)

"""
$(TYPEDSIGNATURES)

In-place partial derivative of `u` along the second dimension. Uses central differences in
the interior; boundary behaviour is controlled by `idx` (default: `FlatIndexing`,
one-sided differences). GPU-compatible and asynchronous (see
[Asynchronous kernel launches](@ref)).
"""
function ∂y!(du, u, dy, idx::AbstractIndexing = FlatIndexing(1, size(u, 2)))
    backend = get_backend(u)
    if backend isa KernelAbstractions.CPU
        kernel! = _∂x₂_col!(backend)
        kernel!(du, u, dy, idx; ndrange = size(u, 2))
    else
        kernel! = _∂y!(backend)
        kernel!(du, u, dy, idx; ndrange = size(u))
    end
    return nothing
end

∂y!(du, u, dy, ::Integer) = ∂y!(du, u, dy)

"""
$(TYPEDSIGNATURES)

In-place partial derivative of `u` along the third dimension. Uses central differences in
the interior; boundary behaviour is controlled by `idx` (default: `FlatIndexing`,
one-sided differences). GPU-compatible and asynchronous (see
[Asynchronous kernel launches](@ref)).
"""
function ∂x₃!(du, u, dz, idx::AbstractIndexing = FlatIndexing(1, size(u, 3)))
    backend = get_backend(u)
    kernel! = _∂x₃!(backend)
    kernel!(du, u, dz, idx; ndrange = size(u))
    return nothing
end

"""
$(TYPEDSIGNATURES)

In-place vertical derivative of `u` in a terrain-following sigma coordinate system.

The effective physical spacing at level `k` is `H[i,j] · (ζ_aa[kp1] - ζ_aa[km1])`,
where `ζ_aa` are the sigma midpoint positions given by `transform` and `H` is the local
ice thickness. Boundary behaviour along the vertical is controlled by `idx`
(default: `FlatIndexing`, one-sided differences). GPU-compatible and asynchronous (see
[Asynchronous kernel launches](@ref)).
"""
function ∂x₃!(du, u, H::AbstractMatrix, transform::AbstractSigmaTransform,
               idx::AbstractIndexing = FlatIndexing(1, size(u, 3)))
    T    = eltype(u)
    ζ_aa = similar(u, size(u, 3))
    copyto!(ζ_aa, get_ζ_aa(T, transform))
    backend = get_backend(u)
    kernel! = _∂x₃_sigma!(backend)
    kernel!(du, u, ζ_aa, H, idx; ndrange = size(u))
    return nothing
end

"""
$(TYPEDSIGNATURES)

Partial derivative of `u` along the first dimension. See [`∂x!`](@ref).
"""
∂x₁(u, dx, idx::AbstractIndexing = FlatIndexing(1, size(u, 1))) =
    (du = similar(u); ∂x!(du, u, dx, idx); du)

"""
$(TYPEDSIGNATURES)

Partial derivative of `u` along the second dimension. See [`∂y!`](@ref).
"""
∂x₂(u, dy, idx::AbstractIndexing = FlatIndexing(1, size(u, 2))) =
    (du = similar(u); ∂y!(du, u, dy, idx); du)

"""
$(TYPEDSIGNATURES)

Partial derivative of `u` along the third dimension. See [`∂x₃!`](@ref).
"""
∂x₃(u, dz, idx::AbstractIndexing = FlatIndexing(1, size(u, 3))) =
    (du = similar(u); ∂x₃!(du, u, dz, idx); du)

"""
$(TYPEDSIGNATURES)

Vertical derivative of `u` in sigma coordinates. See [`∂x₃!`](@ref).
"""
∂x₃(u, H::AbstractMatrix, transform::AbstractSigmaTransform,
    idx::AbstractIndexing = FlatIndexing(1, size(u, 3))) =
    (du = similar(u); ∂x₃!(du, u, H, transform, idx); du)

"""
$(TYPEDSIGNATURES)

In-place computation of both planar partial derivatives of `u`. On GPU backends a single
fused kernel reads `u` only once, saving memory bandwidth. On CPU the two optimised
column kernels ([`∂x!`](@ref), [`∂y!`](@ref)) are called sequentially; they already
achieve good SIMD efficiency individually and the working set fits in cache. Boundary
behaviour is controlled by `idx₁` and `idx₂` independently (default: `FlatIndexing`).
GPU-compatible and asynchronous (see [Asynchronous kernel launches](@ref)).
"""
function ∂x₁₂!(du₁, du₂, u, dx, dy,
    idx₁::AbstractIndexing = FlatIndexing(1, size(u, 1)),
    idx₂::AbstractIndexing = FlatIndexing(1, size(u, 2)))
    backend = get_backend(u)
    if backend isa KernelAbstractions.CPU
        ∂x!(du₁, u, dx, idx₁)
        ∂y!(du₂, u, dy, idx₂)
    else
        kernel! = _∂x₁₂!(backend)
        kernel!(du₁, du₂, u, dx, dy, idx₁, idx₂; ndrange = size(u))
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Fused computation of both planar partial derivatives of `u`. See [`∂x₁₂!`](@ref).
"""
function ∂x₁₂(u, dx, dy,
    idx₁::AbstractIndexing = FlatIndexing(1, size(u, 1)),
    idx₂::AbstractIndexing = FlatIndexing(1, size(u, 2)))
    du₁ = similar(u)
    du₂ = similar(u)
    ∂x₁₂!(du₁, du₂, u, dx, dy, idx₁, idx₂)
    return du₁, du₂
end

@kernel function _∂x!(du, u, dx, idx)
    i, j = @index(Global, NTuple)
    im1, ip1, h = stencil_fd(i, idx)
    @inbounds du[i, j] = (u[ip1, j] - u[im1, j]) / (h * dx)
end

# CPU-optimised variant: dispatches one thread per column (j), then loops over i
# internally.  The stencil_fd call is loop-invariant for the boundary elements,
# and the interior uses a scalar reciprocal so LLVM emits SIMD multiplications
# rather than per-lane vector divisions.
@kernel function _∂x₁_col!(du, u, dx, idx)
    j  = @index(Global)
    i1 = idx.i1
    i2 = idx.i2
    # Lower boundary (one-sided or reflected depending on idx)
    im1, ip1, h = stencil_fd(i1, idx)
    @inbounds du[i1, j] = (u[ip1, j] - u[im1, j]) / (h * dx)
    # Interior: h = 2 always → hoist reciprocal, enable SIMD multiply
    inv_2dx = inv(2 * dx)
    @inbounds for i in (i1 + 1):(i2 - 1)
        du[i, j] = (u[i + 1, j] - u[i - 1, j]) * inv_2dx
    end
    # Upper boundary (skip when domain has only one point)
    if i2 > i1
        im1, ip1, h = stencil_fd(i2, idx)
        @inbounds du[i2, j] = (u[ip1, j] - u[im1, j]) / (h * dx)
    end
end

@kernel function _∂y!(du, u, dy, idx)
    i, j = @index(Global, NTuple)
    jm1, jp1, h = stencil_fd(j, idx)
    @inbounds du[i, j] = (u[i, jp1] - u[i, jm1]) / (h * dy)
end

# CPU-optimised variant: dispatches one thread per column (j) and loops over i
# internally.  stencil_fd(j, idx) is computed once per dispatch (h is uniform
# across all i for a given j), so the entire inner loop uses a scalar reciprocal
# and LLVM emits SIMD multiplications rather than per-lane vector divisions.
# No boundary peeling is needed because h depends only on j (the dispatch index).
@kernel function _∂x₂_col!(du, u, dy, idx)
    j = @index(Global)
    jm1, jp1, h = stencil_fd(j, idx)
    inv_h_dy = inv(h * dy)
    @inbounds for i in axes(u, 1)
        du[i, j] = (u[i, jp1] - u[i, jm1]) * inv_h_dy
    end
end

@kernel function _∂x₃!(du, u, dz, idx)
    i, j, k = @index(Global, NTuple)
    km1, kp1, h = stencil_fd(k, idx)
    @inbounds du[i, j, k] = (u[i, j, kp1] - u[i, j, km1]) / (h * dz)
end

# Sigma-coordinate variant: non-uniform ζ spacing, physical dz = H[i,j] * (ζ_aa[kp1]-ζ_aa[km1]).
@kernel function _∂x₃_sigma!(du, u, ζ_aa, H, idx)
    i, j, k = @index(Global, NTuple)
    km1, kp1, _ = stencil_fd(k, idx)
    @inbounds du[i, j, k] = (u[i, j, kp1] - u[i, j, km1]) / ((ζ_aa[kp1] - ζ_aa[km1]) * H[i, j])
end

@kernel function _∂x₁₂!(du₁, du₂, u, dx, dy, idx₁, idx₂)
    i, j = @index(Global, NTuple)
    im1, ip1, h₁ = stencil_fd(i, idx₁)
    jm1, jp1, h₂ = stencil_fd(j, idx₂)
    @inbounds begin
        du₁[i, j] = (u[ip1, j] - u[im1, j]) / (h₁ * dx)
        du₂[i, j] = (u[i, jp1] - u[i, jm1]) / (h₂ * dy)
    end
end


# @dev TODO: horizontal/time derivatives in a sigma-coordinate system are not yet
# implemented; they would take the form:
#    ∂x₁(u, ∂ζ_∂x₁) = ∂x₁(u)|_σ + ∂ζ_∂x₁ * ∂ζ(u)
#    ∂x₂(u, ∂ζ_∂x₂) = ∂x₂(u)|_σ + ∂ζ_∂x₂ * ∂ζ(u)
#    ∂x₃(u, ∂ζ_∂x₃) = ∂ζ_∂x₃ * ∂ζ(u)
#    ∂t(u,  ∂ζ_∂t)  = ∂t(u)|_σ  + ∂ζ_∂t  * ∂ζ(u)

"""
$(TYPEDSIGNATURES)

An abstract type for indexing strategies used in finite difference stencils. See
[`stencil_fd`](@ref) and the various `index` methods.

# Available subtypes:
 - [`StrictIndexing`](@ref)
 - [`FlatIndexing`](@ref)
 - [`ReflectiveIndexing`](@ref)
 - [`PeriodicIndexing`](@ref)
"""
abstract type AbstractIndexing end

"""
$(TYPEDSIGNATURES)

A struct representing strict indexing, where the caller guarantees that all indices are within bounds. No runtime checks are performed, and the effective stencil width is always 2 (central difference). GPU-compatible.
"""
struct StrictIndexing <: AbstractIndexing
    i1::Int
    i2::Int
end

"""
$(TYPEDSIGNATURES)

A struct representing flat indexing, where out-of-bounds indices are clamped to the nearest valid index. This gives one-sided differences at boundaries and central differences in the interior. GPU-compatible.
"""
struct FlatIndexing <: AbstractIndexing
    i1::Int
    i2::Int
end
"""
$(TYPEDSIGNATURES)

A struct representing reflective indexing, where out-of-bounds indices are reflected back into the domain. This gives zero-gradient boundary conditions. GPU-compatible.
"""
struct ReflectiveIndexing <: AbstractIndexing
    i1::Int
    i2::Int
end
"""
$(TYPEDSIGNATURES)

A struct representing periodic indexing, where out-of-bounds indices wrap around to the opposite side of the domain. GPU-compatible.
"""
struct PeriodicIndexing <: AbstractIndexing
    i1::Int
    i2::Int
end

"""
$(TYPEDSIGNATURES)

Return the index obtained by shifting `i` by `d` according to the indexing strategy `idx`.
"""
@inline index(i, d, ::StrictIndexing) = i + d
# No runtime bounds check: caller is responsible for staying within [i1, i2].
# error() is not GPU-safe, so StrictIndexing simply performs the unchecked shift.

@inline index(i, d, idx::FlatIndexing) = clamp(i + d, idx.i1, idx.i2)
@inline function index(i, d, idx::ReflectiveIndexing)
    j = i + d
    j < idx.i1 && return idx.i1 + (idx.i1 - j)
    j > idx.i2 && return idx.i2 - (j - idx.i2)
    return j
end
@inline function index(i, d, idx::PeriodicIndexing)
    j = i + d
    j < idx.i1 && return idx.i2 - (idx.i1 - j - 1)
    j > idx.i2 && return idx.i1 + (j - idx.i2 - 1)
    return j
end

"""
    stencil_fd(i, idx) -> (im1, ip1, h)

Return the backward index `im1`, forward index `ip1`, and effective stencil width `h`
(in units of one grid spacing) for evaluating `(u[ip1] - u[im1]) / (h * dx)`.

- `FlatIndexing`: `h = ip1 - im1`, giving 1 (one-sided) at boundaries and 2 (central)
  in the interior.
- `ReflectiveIndexing`: `h = 2` always; zero gradient falls out naturally since
  `im1 == ip1` at the boundary.
- `PeriodicIndexing`: `h = 2` always (central difference wrapping around).
- `StrictIndexing`: `h = 2` always; caller guarantees the point is interior.
"""
@inline function stencil_fd(i, idx::FlatIndexing)
    im1 = index(i, -1, idx)
    ip1 = index(i, +1, idx)
    return im1, ip1, ip1 - im1
end
@inline function stencil_fd(i, idx::ReflectiveIndexing)
    im1 = index(i, -1, idx)
    ip1 = index(i, +1, idx)
    return im1, ip1, 2
end
@inline function stencil_fd(i, idx::PeriodicIndexing)
    im1 = index(i, -1, idx)
    ip1 = index(i, +1, idx)
    return im1, ip1, 2
end
@inline function stencil_fd(i, idx::StrictIndexing)
    im1 = index(i, -1, idx)
    ip1 = index(i, +1, idx)
    return im1, ip1, 2
end

"""
$(TYPEDSIGNATURES)

Convenience methods for getting the full stencil of indices needed to compute finite differences in multiple dimensions. For example, `stencil(i, j, idx_i, idx_j)` returns `(im1_i, ip1_i, im1_j, ip1_j)`, which can be used to compute the x and y derivatives at `(i, j)`. GPU-compatible.
"""
function stencil(i, idx)
    return index(i, -1, idx), index(i, 1, idx)
end
function stencil(i, j, i_idx, j_idx)
    return stencil(i, i_idx)..., stencil(j, j_idx)...
end
function stencil(i, j, k, i_idx, j_idx, k_idx)
    return stencil(i, i_idx)..., stencil(j, j_idx)..., stencil(k, k_idx)...
end

"""
$(TYPEDSIGNATURES)

Convert between 2D and 1D indices for x-flattened arrays. Necessary for assembly of linear problem.
"""
function _ij2n_ux(i, j, ::Integer, ny)
    return (i - 1) * ny + j
end

"""
$(TYPEDSIGNATURES)

Convert between 2D and 1D indices for y-flattened arrays. Necessary for assembly of linear problem.
"""
function _ij2n_uy(i, j, nx, ny)
    return (i - 1) * ny + j + nx * ny
end

#=
function von_neumann_neighbours(I::CartesianIndex{N}) where N
    neighbours = CartesianIndex[]
    for d in 1:N
        push!(neighbours, I + CartesianIndex(ntuple(i -> i == d ? -1 : 0, N)))
        push!(neighbours, I + CartesianIndex(ntuple(i -> i == d ? 1 : 0, N)))
    end
    return neighbours
end

function moore_neighbours(I::CartesianIndex{N}) where N
    neighbours = CartesianIndex[]
    for offset in Iterators.product(ntuple(_ -> -1:1, N)...)
        if any(x -> x != 0, offset)
            push!(neighbours, I + CartesianIndex(offset))
        end
    end
    return neighbours
end
=#

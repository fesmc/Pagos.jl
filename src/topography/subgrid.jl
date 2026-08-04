abstract type AbstractGroundedFraction end

struct GridGroundedFraction <: AbstractGroundedFraction end

struct SubgridGroundedFraction <: AbstractGroundedFraction
    H_ice_neighbours::V1        # Vector{Real}
    mask_ice_neighbours::V2     # Vector{Bool}
end

function SubgridGroundedFraction(; dimsize = 2, T = Float32)
    H_ice_neighbours = zeros(T, 2 ^ dimsize)
    mask_ice_neighbours = falses(Bool, 2 ^ dimsize)
    return SubgridGroundedFraction(H_ice_neighbours, mask_ice_neighbours)
end

function ice_fraction(H_ice::T) where {T<:Real}
    if H_ice > 0
        f_ice = 1
    else
        f_ice = 0
    end
    return f_ice
end

function ice_fraction(H_ice, ggf::GridGroundedFraction)
    return ice_fraction.(H_ice)
end


function ice_neighbours!(H_ice_eff, n_ice_neighbours, H_ice, f_ice)
    neighbours = von_neumann_neighbours(CartesianIndex(size(H_ice) .÷ 2))  # dummy index to get size
    H_ice_neighbours = [H_ice[J] for J in neighbours]


    for I in CartesianIndices(H_ice)
        neighbours .= von_neumann_neighbours(I)
        H_ice_neighbours .= [H_ice[J] for J in neighbours]
        n_ice_neighbours[I] = sum(f_ice[J] for J in neighbours)
    end
end

function ice_fraction(H_ice, ggf::SubgridGroundedFraction)
    f_ice = ice_fraction.(H_ice)

    return f_ice
end

"""
    update_grounded_fraction!(fg_ac, Hg, idx)

Update the grounded fraction of the ice shelf as in Eq. (5) of [robinson-description-2020](@citet).
"""
function update_grounded_fraction!(fg_ac, Hg, idx)
    return nothing
end
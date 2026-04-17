
"""
$(TYPEDSIGNATURES)

Compute the effective ice thickness `H_eff`:

```math
\begin{aligned}
H_{eff} = \begin{cases}\dfrac{H_{ice}}{f_{ice}} & \text{for } f_{ice} > 0 \\ 0 & \text{for } f_{ice} = 0 \end{cases}
\end{aligned}
```
"""
function H_ice_effective(H_ice, f_ice)
    if f_ice > 0
        H_eff = H_ice / f_ice
    else
        H_eff = 0
    end
    return H_eff
end

function height_above_floatation(H, ρ_ratio, z_sl, z_bed)
    return H - ρ_ratio * max(z_sl - z_bed, 0)
end


function seawater_depth(ρ_seawater__div__ρ_ice, H_ocn)
    if H_ocn > 0
        H_ocn_now = H_ocn * ρ_seawater__div__ρ_ice
    else
        H_ocn_now = 0.0
    end
    return H_ocn_now
end

struct Densities{T}
    ρ_ice::T
    ρ_seawater::T
    ρ_ice__tim__g::T
    ρ_ice__div__ρ_seawater::T
    ρ_seawater__div__ρ_ice::T
end

function surface_elevation(H_grnd, H_eff, ρ_ratio, z_sl, z_bed)
    if H_grnd > 0
        z_srf = z_bed + H_eff
    else
        z_srf = z_sl + (1 - ρ_ratio) * H_eff
    end
    return z_srf
end

function maximal_surface_elevation()
    return max(z_bed + H_eff, z_sl + (1 - ρ_ratio) * H_eff)
end


function H_ice_grounded(H_eff, z_sl, z_bed, ρ_ratio)
    H_grnd = height_above_floatation(H_eff, ρ_ratio, z_sl, z_bed)
    return H_grnd
end

abstract type AbstractGroundedFraction end

struct GridGroundedFraction <: AbstractGroundedFraction
end

struct SubgridGroundedFraction <: AbstractGroundedFraction
    H_ice_neighbours::V1        # Vector{Real}
    mask_ice_neighbours::V2     # Vector{Bool}
end

function SubgridGroundedFraction(; dimsize = 2, T = Float32)
    H_ice_neighbours = zeros(T, 2 ^ dimsize)
    mask_ice_neighbours = falses(Bool, 2 ^ dimsize)
    return SubgridGroundedFraction(H_ice_neighbours, mask_ice_neighbours)
end

function ice_fraction(H_ice::T) where T<:Real
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
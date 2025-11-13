###########################################################
# Structs
###########################################################

"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the flow law computation via [`viscosity`](@ref) and [`viscosity!`](@ref).
"""
abstract type AbstractFlowLaw end

"""
$(TYPEDSIGNATURES)

Flow law with constant viscosity.
# Fields
 - `η::T`: constant viscosity.
"""
struct ConstantViscosityFlowLaw{T} <: AbstractFlowLaw
    η::T
end

"""
$(TYPEDSIGNATURES)

Flow law combining a rate factor and a creep function.
"""
struct RateCreepFlowLaw{
    RF,     # <: AbstractRateFactor,
    CF,     # <: AbstractCreepFunction,
} <: AbstractFlowLaw
    rate_factor::RF
    creep_function::CF
end

"""
$(TYPEDSIGNATURES)

Convenience function to create a Glen-Nye flow law with [`ArrheniusRateFactor`](@ref) and [`GlenNyeCreepFunction`](@ref).
"""
function GlenNyeFlowLaw()
    rf = ArrheniusRateFactor()
    cf = GlenNyeCreepFunction()
    return RateCreepFlowLaw(rf, cf)
end

"""
$(TYPEDSIGNATURES)

Convenience function to create a Regularized Glen-Nye flow law with [`ArrheniusRateFactor`](@ref) and [`RegularizedGlenNyeCreepFunction`](@ref).
"""
function RegularizedGlenNyeFlowLaw()
    rf = ArrheniusRateFactor()
    cf = RegularizedGlenNyeCreepFunction()
    return RateCreepFlowLaw(rf, cf)
end

"""
$(TYPEDSIGNATURES)

Convenience function to create a Smith-Morland flow law with [`SmithMorlandRateFactor`](@ref) and [`SmithMorlandCreepFunction`](@ref).
"""
function SmithMorlandFlowLaw()
    rf = SmithMorlandRateFactor()
    cf = SmithMorlandCreepFunction()
    return RateCreepFlowLaw(rf, cf)
end

###########################################################
# Dispatch
###########################################################

"""
$(TYPEDSIGNATURES)

Get the viscosity `η` based on the rate factor `A`, creep function `f`, and flow law parameterization `flowlaw<:AbstractFlowLaw`.
"""
function viscosity(
    A,
    f,
    flowlaw::ConstantViscosityFlowLaw,
)
    return flowlaw.η
end

function viscosity(
    A,  # rate factor
    f,  # creep function
    flowlaw::RateCreepFlowLaw,
)
    return 0.5 ./ (A .* f)
end

function viscosity(
    A::M1,
    f::M2,
    flowlaw,
) where {M1<:AbstractArray, M2<:AbstractArray}
    η = similar(A)
    viscosity!(η, A, f, flowlaw)
    return η
end

"""
$(TYPEDSIGNATURES)

Same as [`viscosity`](@ref) but operates in place.
"""
function viscosity!(η, A, f, flowlaw)
    map!((a, s) -> viscosity(a, s, flowlaw), η, A, f)
    return nothing
end

# function update_ice_viscosity!(
#     η::Array{T, 3},
#     A::Array{T, 3},
#     σ_e::Array{T, 3},
#     flowlaw::ConstantViscosityFlowLaw{T},
#     idx::Matrix{CartesianIndex{2}},
# ) where {T<:AbstractFloat}

#     for i in idx
#         for l in axes(η, 3)
#             η[i, l] .= flowlaw.η
#         end
#     end
#     return
# end

# function update_ice_viscosity!(
#     η::Array{T, 3},
#     A::Array{T, 3},
#     σ_e::Array{T, 3},
#     flowlaw::GlenViscosity{T},
#     idx::Matrix{CartesianIndex{2}},
# ) where {T<:AbstractFloat}

#     for i in idx
#         for l in axes(η, 3)
#             η[i, l] = 0.5 / (A[i, l] * σ_e[i, l]^(flowlaw.n-1))
#         end
#     end
# end

# function update_ice_viscosity!(
#     η::Array{T, 3},
#     A::Array{T, 3},
#     σ_e::Array{T, 3},
#     flowlaw::RegularizedGlenViscosity{T},
#     idx::Matrix{CartesianIndex{2}},
# ) where {T<:AbstractFloat}

#     (; n, σ_0) = flowlaw
#     for i in idx
#         for l in axes(η, 3)
#             η[i, l] = 0.5 / (A[i, l] * (σ_e[i, l]^(n-1) + σ_0^(n-1)))
#         end
#     end
# end







function calc_visc_eff_2D_aa(
    ux,
    uy,
    ATT,
    f_ice,
    dx,
    dy;
    n_glen = 3,
    eps_0 = 1e-6,
)
    # Calculate 3D effective viscosity following L19, Eq. 2
    # Use of eps_0 ensures non-zero positive viscosity value everywhere 
    # Note: viscosity is first calculated on ab-nodes, then 
    # unstaggered back to aa-nodes. This ensures more stability for 
    # visc_eff (less likely to blow up for low strain rates). 

    visc_min = 1e5

    nx, ny = size(ux)

    # Calculate exponents 
    p1 = (1.0 - n_glen) / (2.0 * n_glen)
    p2 = -1.0 / n_glen

    # Calculate squared minimum strain rate 
    eps_0_sq = eps_0 * eps_0

    # Calculate visc_eff on aa-nodes
    visc = fill(visc_min, nx, ny)
    eps_aa = fill(eps_0_sq, nx, ny)

    for i = 1:nx
        for j = 1:ny

            if f_ice[i, j] == 1.0

                im1, ip1, jm1, jp1 = periodic_indices(i, j, nx, ny)

                # Get strain rate terms
                dudx_aa = (ux[i, j] - ux[im1, j]) / dx
                dvdy_aa = (uy[i, j] - uy[i, jm1]) / dy

                dudy_aa_1 = (ux[i, jp1] - ux[i, jm1]) / (2.0 * dy)
                dudy_aa_2 = (ux[im1, jp1] - ux[im1, jm1]) / (2.0 * dy)
                dudy_aa = 0.5 * (dudy_aa_1 + dudy_aa_2)

                dvdx_aa_1 = (uy[ip1, j] - uy[im1, j]) / (2.0 * dx)
                dvdx_aa_2 = (uy[ip1, jm1] - uy[im1, jm1]) / (2.0 * dx)
                dvdx_aa = 0.5 * (dvdx_aa_1 + dvdx_aa_2)

                # Calculate the total effective strain rate from L19, Eq. 21 
                eps_sq_aa =
                    dudx_aa^2 +
                    dvdy_aa^2 +
                    dudx_aa * dvdy_aa +
                    0.25 * (dudy_aa + dvdx_aa)^2 +
                    eps_0_sq
                eps_aa[i, j] = sqrt(eps_sq_aa)

                # Get rate factor on central node
                ATT_aa = ATT[i, j]

                # Calculate effective viscosity on ab-nodes
                visc[i, j] = 0.5 * (eps_sq_aa)^(p1) * ATT_aa^(p2)

            end
        end
    end

    #println("eps: ", extrema(eps_aa))

    return visc

end

"""
    calc_visc_eff_2D_nodes(ux,uy,ATT,H_ice,f_ice,dx,dy,xn,yn;n_glen=3,eps_0=1e-6,wtn=fill(1.0,length(xn)))

Calculate 3D effective viscosity following L19, Eq. 2
Use of eps_0 ensures non-zero positive viscosity value everywhere 
Note: viscosity is first calculated on ab-nodes, then 
unstaggered back to aa-nodes. This ensures more stability for 
visc_eff (less likely to blow up for low strain rates). 

Given ux on acx-nodes and uy on acy-nodes, get both quantities 
on node locations of choice [xn;yn]. Viscosity will be calculated
at those locations and the desired weighting wtn will be applied to each node.
"""
function calc_visc_eff_2D_nodes(
    ux,
    uy,
    ATT,
    f_ice,
    dx,
    dy;
    n_glen = 3,
    eps_0 = 1e-6,
)

    visc_min = 1e5
    nx, ny = size(ux)

    # Calculate exponents 
    p1 = (1.0 - n_glen) / (2.0 * n_glen)
    p2 = -1.0 / n_glen

    # Calculate squared minimum strain rate 
    eps_0_sq = eps_0 * eps_0

    # Populate strain rates over the whole domain on acx- and acy-nodes
    dudx = fill(0.0, nx, ny)
    dvdy = fill(0.0, nx, ny)
    dudy = fill(0.0, nx, ny)
    dvdx = fill(0.0, nx, ny)

    for i = 1:nx
        for j = 1:ny
            im1, ip1, jm1, jp1 = periodic_indices(i, j, nx, ny)
            dudx[i, j] = (ux[ip1, j] - ux[im1, j]) / (2.0 * dx)
            dudy[i, j] = (ux[i, jp1] - ux[i, jm1]) / (2.0 * dy)
            dvdx[i, j] = (uy[ip1, j] - uy[im1, j]) / (2.0 * dx)
            dvdy[i, j] = (uy[i, jp1] - uy[i, jm1]) / (2.0 * dy)
        end
    end

    # Calculate visc_eff on aa-nodes

    visc = fill(visc_min, nx, ny)
    eps_aa = fill(eps_0_sq, nx, ny)

    wt0 = 1.0 / sqrt(3)
    xn = [wt0, -wt0, -wt0, wt0]
    yn = [wt0, wt0, -wt0, -wt0]
    wtn = [1.0, 1.0, 1.0, 1.0]

    for i = 1:nx
        for j = 1:ny

            if f_ice[i, j] == 1.0

                im1, ip1, jm1, jp1 = periodic_indices(i, j, nx, ny)

                # Get strain rate terms on node locations
                dudxn = acx_to_nodes(dudx, i, j, xn, yn)
                dudyn = acx_to_nodes(dudy, i, j, xn, yn)

                dvdxn = acy_to_nodes(dvdx, i, j, xn, yn)
                dvdyn = acy_to_nodes(dvdy, i, j, xn, yn)

                # Calculate the total effective strain rate from L19, Eq. 21 
                eps_sq_n =
                    dudxn .^ 2 + dvdyn .^ 2 .+ dudxn .* dvdyn .+
                    0.25 .* (dudyn .+ dvdxn) .^ 2 .+ eps_0_sq
                eps_aa = sum(sqrt.(eps_sq_n)) / length(eps_sq_n)

                # Get rate factor on central node
                ATT_aa = ATT[i, j]

                # Calculate effective viscosity on ab-nodes
                viscn = 0.5 .* (eps_sq_n) .^ (p1) .* ATT_aa^(p2)
                visc[i, j] = sum(viscn .* wtn) / sum(wtn)

            end
        end
    end
    #println("eps: ", extrema(eps_aa))
    return visc
end
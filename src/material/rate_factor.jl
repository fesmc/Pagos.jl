###########################################################
# Structs
###########################################################

"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the rate factor computation via [`rate_factor`](@ref).
"""
abstract type AbstractRateFactor{T<:AbstractFloat} end

"""
$(TYPEDSIGNATURES)

Rate factor for ice viscosity following a constant value.

# Fields
- `A::T=1e-16`: rate factor.
"""
@kwdef struct ConstantRateFactor{T} <: AbstractRateFactor{T}
    A::T = 1e-16
end

"""
$(TYPEDSIGNATURES)

Rate factor for ice viscosity following Arrhenius' law.

# Fields
 - `E_f::T=1.0`: enhancement factor.
 - `T_p1_p2::T=263.15`: breakpoint temperature following Patterson (1994).
 - `A_0_p1::T=3.985e-13`: "s^-1 Pa^-3" pre-exponential factor.
 - `A_0_p2::T=1.916e3`: piecewise definition following Patterson (1994).
 - `Q_a_p1::T=60e3`: "J mol^-1" activation energy.
 - `Q_a_p2::T=139e3`: piecewise definition following Patterson (1994).
"""
@kwdef struct ArrheniusRateFactor{T} <: AbstractRateFactor{T}
    E_f::T = 1.0
    T_p1_p2::T = 263.15     # breakpoint temperature following Patterson (1994)
    A_0_p1::T = 3.985e-13   # "s^-1 Pa^-3" pre-exponential factor
    A_0_p2::T = 1.916e3     # piecewise definition following Patterson (1994)
    Q_a_p1::T = 60e3        # "J mol^-1" activation energy
    Q_a_p2::T = 139e3       #  piecewise definition following Patterson (1994)
    R::T = 8.314            # "J K^-1 mol^-1" universal gas constant
end

"""
$(TYPEDSIGNATURES)

Warning: This is not fully supported yet, since the input fields are different than those used in Glen's flow law, which is much more common.

Rate factor for ice viscosity following Smith and Morland (1981).

# Fields
 - `p1::T=0.7242`
 - `p2::T=0.3438`
 - `e1::T=11.9567`
 - `e2::T=2.9494`
"""
@kwdef struct SmithMorlandRateFactor{T} <: AbstractRateFactor{T}
    p1::T = 0.7242
    p2::T = 0.3438
    e1::T = 11.9567
    e2::T = 2.9494
    T_0::T = 273.15
    ΔT::T = 20.0
end

###########################################################
# Functions
###########################################################

"""
$(TYPEDSIGNATURES)

Get the rate factor `A` based on the temperature relative to the pressure melt point `T_relative` and the rate factor parameterization `arf<:AbstractRateFactor`.
"""
function rate_factor(
    T_relative::T,
    arf::ArrheniusRateFactor,
) where {T<:Real}

    (; E_f, T_p1_p2, A_0_p1, A_0_p2, Q_a_p1, Q_a_p2, R) = arf
    if T_relative <= T_p1_p2
        A = E_f * A_0_p1 * exp(-Q_a_p1 / (R * T_relative))
    else
        A = E_f * A_0_p2 * exp(-Q_a_p2 / (R * T_relative))
    end
    return A
end

function rate_factor(
    T_relative::T,
    smr::SmithMorlandRateFactor,
) where {T<:Real}

    (; p1, p2, e1, e2) = smr
    T_bar = (T - T_0) / ΔT
    A = p1 * exp(e1 * T_relative) + p2 * exp(e2 * T_relative)
    return A
end

function rate_factor(
    T_relative::M,
    arf::ARF,
) where {M<:AbstractArray, ARF<:AbstractRateFactor}
    A = similar(T_relative)
    rate_factor!(A, T_relative, arf)
    return A
end

"""
$(TYPEDSIGNATURES)

Get the rate factor `A` based on the temperature relative to the pressure melt point `T_relative` and the rate factor parameterization `arf<:AbstractRateFactor`. Optionally, a mask can be provided to only compute the rate factor for specific indices.
"""
function rate_factor!(
    A,
    T_relative,
    arf::AbstractRateFactor,
)
    map!(x -> rate_factor(x, arf), A, T_relative)
    return
end

function rate_factor!(
    A,
    T_relative,
    arf::AbstractRateFactor,
    mask,
)
    map!( x -> rate_factor(x, arf), view(A, mask), view(T_relative, mask))
    return
end


# """
#     update_rate_factor!(A, T_relative, ae, idx)

# Update the rate factor `A` based on the temperature relative to the pressure melt
# point `T_relative` and the rate factor `ae<:AbstractRateFactor`.
# """
# function update_rate_factor!(
#     A::Array{T, 3},
#     T_relative::Array{T, 3},
#     cae::ConstantRateFactor{T},
#     idx::Matrix{CartesianIndex{2}},
# ) where {T<:AbstractFloat}

#     for i in idx
#         for l in axes(A, 3)
#             A[i, l] = cae.A
#         end
#     end
#     return
# end

# function update_rate_factor!(
#     A::Array{T, 3},
#     T_relative::Array{T, 3},
#     aae::ArrheniusRateFactor{T},
#     idx::Matrix{CartesianIndex{2}},
# ) where {T<:AbstractFloat}

#     (; E_f, T_p1_p2, A_0_p1, A_0_p2, Q_a_p1, Q_a_p2, R) = aae
#     for i in idx
#         for l in axes(A, 3)
#             if T_relative[i, l] <= T_p1_p2
#                 A[i, l] = E_f * A_0_p1 * exp(-Q_a_p1 / (R * T_relative[i, l]))
#             else
#                 A[i, l] = E_f * A_0_p2 * exp(-Q_a_p2 / (R * T_relative[i, l]))
#             end
#         end
#     end
# end

# function update_rate_factor!(
#     A::Array{T, 3},
#     T_relative::Array{T, 3},
#     smae::SmithMorlandRateFactor{T},
#     idx::Matrix{CartesianIndex{2}},
# ) where {T<:AbstractFloat}

#     (; p1, p2, e1, e2) = smae
#     for i in idx
#         for l in axes(A, 3)
#             A[i, l] = p1 * exp(e1 * T_relative[i, l]) + p2 * exp(e2 * T_relative[i, l])
#         end
#     end
#     return nothing
# end

"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the sigma transform.

# Available subtypes:
 - [`PowerSigmaTransform`](@ref)
 - [`ArctanSigmaTransform`](@ref)
"""
abstract type AbstractSigmaTransform end

"""
$(TYPEDSIGNATURES)

A sigma transform that uses a power distribution of sigma levels.
"""
struct PowerSigmaTransform{T} <: AbstractSigmaTransform
    n::Int
    exponent::T
end

LinearSigmaTransform(T, n) = PowerSigmaTransform{T}(n, 1)
QuadraticSigmaTransform(T, n) = PowerSigmaTransform{T}(n, 2)

struct ArctanSigmaTransform{T} <: AbstractSigmaTransform{T}
    n::Int
    stretch_factor::T
end

function get_ζ_aa(T, transform::PowerSigmaTransform)
    (; n, exponent) = transform
    ζ_aa = range(0.0, stop = 1.0, length = n) .^ exponent
    return T.(ζ_aa)
end

function get_ζ_ac(ζ_aa)
    n = length(ζ_aa)
    ζ_ac = zeros(eltype(ζ_aa), n + 1)
    for i in 2:n
        ζ_ac[i] = 0.5 * (ζ_aa[i - 1] + ζ_aa[i])
    end
    ζ_ac[n+1] = 1
    return ζ_ac
end

struct VerticalLayering{T, S}
    transform::S
    ζ_aa::Vector{T}
    ζ_ac::Vector{T}
end

function VerticalLayering(T, transform)
    ζ_aa = get_ζ_aa(T, transform)
    ζ_ac = get_ζ_ac(ζ_aa)
    return VerticalLayering(transform, ζ_aa, ζ_ac)
end

struct CorrectedVerticalLayering{T}
    n::Int
    ζ_aa::Vector{T}
    ζ_ac::Vector{T}
end

function CorrectedVerticalLayering(T, transform)
    ζ_ac = range(0.0, stop = 1.0, length = transform.n + 1) .^ transform.exponent
    ζ_aa = zeros(T, transform.n)
    for i in 1:transform.n
        ζ_aa[i] = 0.5 * (ζ_ac[i] + ζ_ac[i + 1])
    end
    return CorrectedVerticalLayering(transform.n, ζ_aa, T.(ζ_ac))
end

sigma(z, b, H) = (z - b) / H

"""
$(TYPEDSIGNATURES)
"""
# function sigma_transform(transform::LinearSigmaTransform)
#     return range(0.0, 1.0; length=transform.n)
# end


# function sigma_transform(transform::ExponentialSigmaTransform)
#     dz = 1 / transform.n
#     return (exp.(dz:dz:1) .- exp(0)) ./ (exp(1) - exp(0))
# end
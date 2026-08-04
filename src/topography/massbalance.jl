abstract type AbstractMassBalance end

abstract type AbstractGroundingZoneMelt end

struct PrescribedGroundingZoneMelt{T} <: AbstractGroundingZoneMelt
    melt_rate::T   # <: AbstractArray or Real
end

function NoGroundingZoneMelt(; T = Float32)
    return PrescribedGroundingZoneMelt{T}(0)
end

# following leguy
struct PartialGroundingZoneMelt{T} <: AbstractGroundingZoneMelt end

# following leguy
struct FullGroundingZoneMelt{T} <: AbstractGroundingZoneMelt end

# following Juarez-Martinez
struct LinearTidePartialGroundingZoneMelt{T} <: AbstractGroundingZoneMelt end

# TODO: following own implementation
struct NonlinearTidePartialGroundingZoneMelt{T} <: AbstractGroundingZoneMelt end

# TODO: try to include topography, hydrology, bed nature, porosity and velocity + something else?
struct MultilinearTidePartialGroundingZoneMelt{T} <: AbstractGroundingZoneMelt end

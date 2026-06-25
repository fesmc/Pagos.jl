abstract type AbstractDamage end

"""
$(TYPEDSIGNATURES)

Prescribe the time-constant ice damage.
"""
struct PrescribedDamage{M} <: AbstractDamage
    damage::M   # <: AbstractArray or Real
end

function NoDamage()
    return PrescribedDamage(0)
end
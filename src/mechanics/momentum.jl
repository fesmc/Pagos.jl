###############################################################
# MomentumBalance
##############################################################

"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the dynamics type via [`velocity`](@ref).

# Available subtypes:
"""
abstract type AbstractMomentumBalance end

"""
$(TYPEDSIGNATURES)

Enforce the ice dynamics to be zero (no flow).
"""
struct NoMomentumBalance <: AbstractMomentumBalance end

struct SIAMomentumBalance <: AbstractMomentumBalance end
struct SSAMomentumBalance <: AbstractMomentumBalance end
struct SIASSAMomentumBalance <: AbstractMomentumBalance end
struct InertialSIASSAMomentumBalance <: AbstractMomentumBalance end
struct DIVAMomentumBalance <: AbstractMomentumBalance end
struct InertialDIVAMomentumBalance <: AbstractMomentumBalance end
struct BlatterPattynMomentumBalance <: AbstractMomentumBalance end
struct StokesMomentumBalance <: AbstractMomentumBalance end
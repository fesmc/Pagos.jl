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

"""
$(TYPEDSIGNATURES)

Shallow Ice Approximation (SIA): the momentum balance is reduced to local vertical
shearing, appropriate for slow, grounded ice far from the margins.
"""
struct SIAMomentumBalance <: AbstractMomentumBalance end

"""
$(TYPEDSIGNATURES)

Shallow Shelf Approximation (SSA): the momentum balance is reduced to depth-independent
membrane stresses, appropriate for floating ice shelves and fast-sliding streams.
"""
struct SSAMomentumBalance <: AbstractMomentumBalance end

"""
$(TYPEDSIGNATURES)

Hybrid SIA+SSA momentum balance, superposing the vertical shearing of the
[`SIAMomentumBalance`](@ref) with the membrane stresses of the [`SSAMomentumBalance`](@ref).
"""
struct SIASSAMomentumBalance <: AbstractMomentumBalance end

"""
$(TYPEDSIGNATURES)

Inertial variant of the [`SIASSAMomentumBalance`](@ref), retaining the acceleration term so
that the velocity is advanced in time rather than solved for at (quasi-)equilibrium.
"""
struct InertialSIASSAMomentumBalance <: AbstractMomentumBalance end

"""
$(TYPEDSIGNATURES)

Depth-Integrated Viscosity Approximation (DIVA): a higher-order momentum balance that
retains vertical shearing through a depth-integrated effective viscosity.
"""
struct DIVAMomentumBalance <: AbstractMomentumBalance end

"""
$(TYPEDSIGNATURES)

Inertial variant of the [`DIVAMomentumBalance`](@ref), retaining the acceleration term so
that the velocity is advanced in time rather than solved for at (quasi-)equilibrium.
"""
struct InertialDIVAMomentumBalance <: AbstractMomentumBalance end

"""
$(TYPEDSIGNATURES)

Blatter-Pattyn higher-order momentum balance, resolving horizontal velocities over the full
column while neglecting the vertical resistive stresses of the full Stokes system.
"""
struct BlatterPattynMomentumBalance <: AbstractMomentumBalance end

"""
$(TYPEDSIGNATURES)

Full Stokes momentum balance, resolving all stress components without the shallowness or
higher-order approximations of the other [`AbstractMomentumBalance`](@ref) subtypes.
"""
struct StokesMomentumBalance <: AbstractMomentumBalance end
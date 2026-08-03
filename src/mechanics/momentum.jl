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

###############################################################
# Dimensionality groupings
###############################################################
#
# The criterion is **the dimensionality of the unknown the momentum solve iterates**, not
# whether the balance has any vertical structure at all. That is what makes the split
# useful for dispatch: it is exactly the question "does this solve live on `grid2d` or on
# `grid`?", which decides node classes, launcher and state fields throughout.
#
# It puts DIVA in the *2D* group, which is the right answer even though DIVA is a
# higher-order balance with a genuine vertical profile: the depth-integrated balance
# (Robinson et al. 2022, Eq. 14) is solved for the depth-averaged `ū`, `v̄`, and `u(z)` is
# reconstructed afterwards from Eq. 16 as a post-solve diagnostic.
#
# The stub balances (`SIA`, `SIASSA`, the `Inertial*` variants, `NoMomentumBalance`) are
# deliberately in neither group until they are ported — an unported balance should hit a
# `MethodError`, not silently inherit a dispatch that was never written for it.

"""
    MomentumBalance2D

The momentum balances whose solve iterates a **depth-integrated** unknown, on `grid2d`:
[`SSAMomentumBalance`](@ref) and [`DIVAMomentumBalance`](@ref).

DIVA belongs here despite resolving vertical shear: its depth-integrated balance is solved
for `ū`/`v̄`, and the 3D profile `u(z)` is reconstructed afterwards rather than iterated.
Both therefore share the same C-grid layout (`acx`/`acy` velocity faces on `grid2d`) and the
same membrane-stress assembly, which is why they share so many method bodies.
"""
const MomentumBalance2D = Union{SSAMomentumBalance, DIVAMomentumBalance}

"""
    MomentumBalance3D

The momentum balances whose solve iterates a genuine **3D** velocity field, on the column
grid: [`BlatterPattynMomentumBalance`](@ref) and [`StokesMomentumBalance`](@ref).

These are also exactly the balances that resolve the vertical velocity `w`, so `∂w/∂z`
(`velocity.z_dz`) is available and `ε̇_zz` is taken from it directly rather than
reconstructed from incompressibility — the property [`raw_strainrate!`](@ref) dispatches on.
"""
const MomentumBalance3D = Union{BlatterPattynMomentumBalance, StokesMomentumBalance}
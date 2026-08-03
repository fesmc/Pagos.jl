"""
$(TYPEDSIGNATURES)

Compute the strain rate tensor in place.
"""
function strainrate!(strainrate, velocity, material, topo, momentum::AbstractMomentumBalance)
    backend = get_backend(strainrate.xx)
    kernel! = _strainrate_kernel!(backend)
    kernel!(strainrate, velocity, material, topo, momentum; ndrange = length(strainrate.xx))
    return nothing
end

# Per-element fallback: triggers when a momentum balance has no specialized method below,
# i.e. the real extension point. Guarding the launcher instead would never catch this.
function strainrate!(strainrate, velocity, material, topo, momentum::AbstractMomentumBalance, I)
    throw(ArgumentError("Unsupported momentum balance type: $(typeof(momentum))"))
end

function strainrate!(strainrate, velocity, material, topo, momentum::SIAMomentumBalance, I)
    # For SIA: depth-averaged viscosity, thickness, strainrate.xz/yz, velocity.*_bar_dz are all 2D
    strainrate.xz[I] = material.viscosity_depthaveraged[I] * topo.thickness[I] * velocity.depthaverage_x_dz[I]
    strainrate.yz[I] = material.viscosity_depthaveraged[I] * topo.thickness[I] * velocity.depthaverage_y_dz[I]
end

function strainrate!(strainrate, velocity, material, topo, momentum::MomentumBalance2D, I)
    # For SSA/DIVA: depth-averaged viscosity, thickness, strainrate.xx/xy/yy, velocity.*_d* are all 2D
    strainrate.xx[I] = 2 * material.viscosity_depthaveraged[I] * topo.thickness[I] * (2 * velocity.x_dx[I] + velocity.y_dy[I])
    strainrate.xy[I] = material.viscosity_depthaveraged[I] * topo.thickness[I] * (velocity.x_dy[I] + velocity.y_dx[I])
    strainrate.yy[I] = 2 * material.viscosity_depthaveraged[I] * topo.thickness[I] * (velocity.x_dx[I] + 2 * velocity.y_dy[I])
end

function strainrate!(strainrate, velocity, material, topo, momentum::BlatterPattynMomentumBalance, I)
    # For Blatter-Pattyn: 3D viscosity, strainrate.xx/xy/yy/xz/yz, velocity.*_d* are all 3D
    strainrate.xx[I] = 2 * material.viscosity[I] * (2 * velocity.x_dx[I] + velocity.y_dy[I])
    strainrate.xy[I] = material.viscosity[I] * (velocity.x_dy[I] + velocity.y_dx[I])
    strainrate.yy[I] = 2 * material.viscosity[I] * (velocity.x_dx[I] + 2 * velocity.y_dy[I])
    strainrate.xz[I] = material.viscosity[I] * velocity.x_dz[I]
    strainrate.yz[I] = material.viscosity[I] * velocity.y_dz[I]
end

@kernel function _strainrate_kernel!(strainrate, velocity, material, topo, momentum::AbstractMomentumBalance)
    I = @index(Global, Linear)
    @inbounds begin
        strainrate!(strainrate, velocity, material, topo, momentum, I)
        strainrate_effective!(strainrate, velocity, momentum, I) 
    end
end

"""
$(TYPEDSIGNATURES)

Compute the effective strain rate in place.
"""
# Per-element fallback: triggers when a momentum balance has no specialized method below.
function strainrate_effective!(strainrate, velocity, momentum::AbstractMomentumBalance, I)
    throw(ArgumentError("Unsupported momentum balance type: $(typeof(momentum))"))
end

function strainrate_effective!(strainrate, velocity, momentum::SIAMomentumBalance, I)
    strainrate.effective[I] = sqrt(1 / 4 * (velocity.depthaverage_x_dz[I] + velocity.depthaverage_y_dz[I]) ^ 2)
end

function strainrate_effective!(strainrate, velocity, momentum::SSAMomentumBalance, I)
    strainrate.effective[I] = sqrt(
        velocity.x_dx[I]^2 + velocity.y_dy[I]^2 +
        velocity.x_dx[I] * velocity.y_dy[I] +
        1 / 4 * (velocity.x_dy[I] + velocity.y_dx[I]) ^ 2
    )
end

# Deliberately *not* [`MomentumBalance2D`](@ref)/[`MomentumBalance3D`](@ref): this union
# cuts across both. The criterion here is "does ε̇_e include the vertical-shear terms
# ¼(u_z² + v_z²)", which DIVA (a 2D-unknown balance, Robinson et al. 2022 Eq. 13) and
# Blatter-Pattyn (a 3D-unknown one) share while SSA does not. Spelling it out keeps that
# distinct grouping visible rather than hiding it behind a name that would suggest
# dimensionality.
function strainrate_effective!(strainrate, velocity, momentum::MB, I) where MB<:Union{DIVAMomentumBalance, BlatterPattynMomentumBalance}
    strainrate.effective[I] = sqrt(
        velocity.x_dx[I]^2 + velocity.y_dy[I]^2 +
        velocity.x_dx[I] * velocity.y_dy[I] +
        1 / 4 * (velocity.x_dy[I] + velocity.y_dx[I]) ^ 2 +
        1 / 4 * velocity.x_dz[I]^2 +
        1 / 4 * velocity.y_dz[I]^2
    )
end


"""
$(TYPEDSIGNATURES)

Compute the **raw** (unscaled) strain-rate tensor
"""
function raw_strainrate!(strainrate, velocity, momentum::AbstractMomentumBalance)
    backend = get_backend(strainrate.xx)
    kernel! = _raw_strainrate_kernel!(backend)
    kernel!(strainrate, velocity, momentum; ndrange = length(strainrate.xx))
    return nothing
end

@kernel function _raw_strainrate_kernel!(strainrate, velocity, momentum::AbstractMomentumBalance)
    I = @index(Global, Linear)
    @inbounds begin
        raw_strainrate!(strainrate, velocity, momentum, I)
        raw_strainrate_effective!(strainrate, momentum, I)
    end
end

# Default (plane / depth-integrated): ε̇_zz reconstructed from incompressibility.
#
# `I` is annotated `::Integer` (it is the flat `@index(Global, Linear)` of the kernel above,
# so this is behaviour-neutral) purely to keep this 4-argument per-element method from being
# ambiguous with the 4-argument staggered `raw_strainrate!(sr, vel, momentum, rt::Runtime)`
# further down. Left untyped, `(Any, Any, MomentumBalance3D, Any)` and
# `(StrainRateState, VelocityState, AbstractMomentumBalance, Runtime)` match the same call
# with neither more specific, and `Runtime` has an empty type intersection with `Integer`,
# so annotating removes the ambiguity without narrowing anything real.
function raw_strainrate!(strainrate, velocity, momentum::AbstractMomentumBalance, I::Integer)
    dxx = velocity.x_dx[I]
    dyy = velocity.y_dy[I]
    strainrate.xx[I] = dxx
    strainrate.yy[I] = dyy
    strainrate.zz[I] = -(dxx + dyy)                            # incompressibility (continuity)
    strainrate.xy[I] = (velocity.x_dy[I] + velocity.y_dx[I]) / 2
    strainrate.xz[I] = (velocity.x_dz[I] + velocity.z_dx[I]) / 2
    strainrate.yz[I] = (velocity.y_dz[I] + velocity.z_dy[I]) / 2
end

# Full-column: ε̇_zz read directly from the resolved vertical velocity gradient.
function raw_strainrate!(strainrate, velocity, momentum::MomentumBalance3D,
                         I::Integer)
    strainrate.xx[I] = velocity.x_dx[I]
    strainrate.yy[I] = velocity.y_dy[I]
    strainrate.zz[I] = velocity.z_dz[I]                        # ∂w/∂z directly
    strainrate.xy[I] = (velocity.x_dy[I] + velocity.y_dx[I]) / 2
    strainrate.xz[I] = (velocity.x_dz[I] + velocity.z_dx[I]) / 2
    strainrate.yz[I] = (velocity.y_dz[I] + velocity.z_dy[I]) / 2
end

"""
$(TYPEDSIGNATURES)

Compute the effective (second-invariant) raw strain rate from the tensor components
already written by [`raw_strainrate!`](@ref) in the same kernel pass.
"""
function raw_strainrate_effective!(strainrate, momentum::AbstractMomentumBalance, I)
    strainrate.effective[I] = sqrt(
        (strainrate.xx[I]^2 + strainrate.yy[I]^2 + strainrate.zz[I]^2) / 2 +
        strainrate.xy[I]^2 + strainrate.xz[I]^2 + strainrate.yz[I]^2
    )
end


"""
$(TYPEDSIGNATURES)

Compute the velocity gradients in x (`v_x_dx, v_y_x`) and y-direction (`v_x_dy, v_y_y`).
The gradients are computed using the central difference scheme. The input velocities
`v_x` and `v_y` are defined on a staggered grid with dimensions `nx` and `ny`.
The grid spacing in x and y-direction is given by `dx` and `dy`.
"""
function velocitygradients!(velocity::VelocityState, dx, dy)
    (; x, y) = velocity
    return velocitygradients!(velocity.x_dx, velocity.x_dy, velocity.y_dx, velocity.y_dy, x, y, dx, dy)
end

function velocitygradients!(v_x_dx, v_x_dy, v_y_dx, v_y_dy, v_dx, v_dy, dx, dy)
    ∂x!(v_x_dx, v_dx, dx)
    ∂y!(v_x_dy, v_dx, dy)
    ∂x!(v_y_dx, v_dy, dx)
    ∂y!(v_y_dy, v_dy, dy)
    return nothing
end

###############################################################
# Chmy-native, C-grid staggered velocity gradients and strain rates
###############################################################
#
# This is where the C-grid layout is supposed to pay for itself, and the check that it
# does is that **no interpolation appears anywhere in the tensor**. Assigning `u` to `acx`,
# `v` to `acy` and `w` to `aa`/z-`Vertex` fixes every gradient by operator algebra:
#
#   ∂u/∂x : acx    → aa       ∂v/∂x : acy    → ab       ∂w/∂x : aa_ac → acx_ac
#   ∂u/∂y : acx    → ab       ∂v/∂y : acy    → aa       ∂w/∂y : aa_ac → acy_ac
#   ∂u/∂z : acx    → acx_ac   ∂v/∂z : acy    → acy_ac   ∂w/∂z : aa_ac → aa
#
# so each strain-rate component is a sum of terms that already live on the *same* node:
# ε̇_xx = ∂u/∂x at `aa`; ε̇_xy = (∂u/∂y + ∂v/∂x)/2 with both at `ab`; ε̇_xz =
# (∂u/∂z + ∂w/∂x)/2 with both at `acx_ac`. Nothing is staggered mid-formula. The
# interpolations that remain are exactly the two places where physics genuinely mixes node
# classes: the viscosity, which lives at `aa` but has to scale off-diagonal strain rates
# elsewhere, and the second invariants, which are cell-centred quantities built from
# off-diagonal components that are not.
#
# Two structural differences from the collocated kernels above, both forced:
#
#  1. **No shared flat `@index(Global, Linear)`.** The components have four different
#     shapes on the C-grid (`aa`, `ab`, `acx_ac`, `acy_ac`), so one linear index addresses
#     a different physical point in each. The kernels below index every field at *its own*
#     `(i, j, k)`, which is the same physical location precisely because each field's index
#     space is anchored to its own node class. The `Launcher`'s `size(grid, Center()) .+ 2`
#     sweep covers every one of those index spaces at once (`aa` needs `1:nx`, `ab` and
#     `acx_ac` need `1:nx+1`, and the sweep provides `0:nx+1`), so a single fused kernel is
#     still possible — it is the flat index that has to go, not the fusion.
#  2. **The invariants need their own launch.** `ε̇_e` at `aa` reads `ε̇_xy` at neighbouring
#     `ab` nodes, so it cannot be computed in the same pass that writes them — the
#     neighbour may not exist yet. The collocated version got away with one pass because
#     every component was at the same point. Hence `raw_strainrate!` then
#     `raw_strainrate_effective!`, in that order.

# Physical vertical derivative on the terrain-following sigma axis: `∂/∂z = (1/H) ∂/∂ζ`,
# matching the legacy `∂x₃!(du, u, H, transform)` convention (`Δu / (Δζ · H)`). `∂z_σ`
# supplies `∂/∂ζ` — Chmy's own `∂z` is wrong on a non-uniform axis, see
# `src/api/sigma_operators.jl`. Unlike the legacy kernel this returns zero rather than
# `Inf`/`NaN` where there is no ice; an `Inf` here propagates into the whole tensor.
@inline _dz_over_H(dζ, H) = H > zero(H) ? dζ / H : zero(dζ)

@kernel inbounds = true function _velocity_gradients!(velocity, H, mask, grid, grid2d, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    u, v, w = velocity.x, velocity.y, velocity.z
    Z = zero(eltype(velocity.x_dx))

    # `node_active` consults only the horizontal part of the node class, so the `_AC`
    # (z-`Vertex`) variants share their activity with `NODE_ACX`/`NODE_ACY`.
    act_aa  = node_active(mask, NODE_AA, i, j)
    act_ab  = node_active(mask, NODE_AB, i, j)
    act_acx = node_active(mask, NODE_ACX, i, j)
    act_acy = node_active(mask, NODE_ACY, i, j)

    velocity.x_dx[I...] = act_aa ? ∂x(u, grid, I...) : Z    # acx   → aa
    velocity.y_dy[I...] = act_aa ? ∂y(v, grid, I...) : Z    # acy   → aa
    velocity.x_dy[I...] = act_ab ? ∂y(u, grid, I...) : Z    # acx   → ab
    velocity.y_dx[I...] = act_ab ? ∂x(v, grid, I...) : Z    # acy   → ab
    velocity.z_dx[I...] = act_acx ? ∂x(w, grid, I...) : Z   # aa_ac → acx_ac
    velocity.z_dy[I...] = act_acy ? ∂y(w, grid, I...) : Z   # aa_ac → acy_ac

    # H is a `grid2d` field, so it is indexed at k = 1 regardless of the sweep's k, and
    # staggered onto the same horizontal face as the derivative it scales.
    velocity.x_dz[I...] = act_acx ?
        _dz_over_H(∂z_σ(u, grid, I...), lerp(H, NODE_ACX, grid2d, i, j, 1)) : Z
    velocity.y_dz[I...] = act_acy ?
        _dz_over_H(∂z_σ(v, grid, I...), lerp(H, NODE_ACY, grid2d, i, j, 1)) : Z
    velocity.z_dz[I...] = act_aa ?
        _dz_over_H(∂z_σ(w, grid, I...), H[i, j, 1]) : Z
end

"""
$(TYPEDSIGNATURES)

Chmy-native [`velocitygradients!`](@ref): fill all nine velocity-gradient fields of
`velocity` from its `x`/`y`/`z` components, each at the node class the corresponding
operator maps to (see the table in the source). `H` is the cell-centred ice thickness on
`rt.grid2d`, needed for the sigma-coordinate vertical scaling `∂/∂z = (1/H) ∂/∂ζ`.

Distinguished from the collocated method by taking a [`Runtime`](@ref) instead of `dx`,
`dy`. Launched on the column grid, since the gradients are column fields.

!!! note "Vertical derivatives use `∂z_σ`, not Chmy's `∂z`"
    On the sigma `Chmy.FunctionAxis` Chmy's own `∂z` scales by the wrong spacing (see
    [`∂z_σ`](@ref)). The horizontal axes are uniform, so `∂x`/`∂y` are unaffected.

!!! warning "No sigma correction on the horizontal derivatives"
    `∂u/∂x` here is taken at constant ζ, not at constant z — the terrain-following
    correction `-(∂z/∂x)/(∂z/∂ζ) · ∂u/∂ζ` is *not* applied. This matches the collocated
    implementation it replaces (which uses plain `∂x!`), so the two are comparable, but it
    is an approximation both share: it is accurate where the surface and bed slopes are
    small, which is the shallow-ice regime these balances assume anyway.
"""
function velocitygradients!(velocity::VelocityState, H, rt::Runtime,
                            mask::AbstractIceMask = NoMask())
    rt.launch(rt.arch, rt.grid,
              _velocity_gradients! => (velocity, H, mask, rt.grid, rt.grid2d))
    return nothing
end

# The depth-averaged counterpart of `_velocity_gradients!`, for the four horizontal
# gradients of `ū`/`v̄` that the SSA/DIVA membrane stress is built from. Same node algebra
# as the column kernel — `∂x` of an `acx` field lands on `aa`, `∂y` of it on `ab` — but
# every field involved is depth-integrated, so this runs on `grid2d` and there is no
# thickness argument: the sigma scaling `∂/∂z = (1/H)∂/∂ζ` that `H` exists for in the
# column kernel has no counterpart here.
@kernel inbounds = true function _depthaverage_velocity_gradients!(velocity, mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    u, v = velocity.depthaverage_x, velocity.depthaverage_y
    Z = zero(eltype(velocity.depthaverage_x_dx))

    act_aa = node_active(mask, NODE_AA, i, j)
    act_ab = node_active(mask, NODE_AB, i, j)

    velocity.depthaverage_x_dx[I...] = act_aa ? ∂x(u, grid, I...) : Z   # acx → aa
    velocity.depthaverage_y_dy[I...] = act_aa ? ∂y(v, grid, I...) : Z   # acy → aa
    velocity.depthaverage_x_dy[I...] = act_ab ? ∂y(u, grid, I...) : Z   # acx → ab
    velocity.depthaverage_y_dx[I...] = act_ab ? ∂x(v, grid, I...) : Z   # acy → ab
end

"""
$(TYPEDSIGNATURES)

Fill the four horizontal gradients of the depth-averaged velocity —
`velocity.depthaverage_x_dx`/`y_dy` (at `aa`) and `depthaverage_x_dy`/`y_dx` (at `ab`) —
from `velocity.depthaverage_x`/`depthaverage_y`. These are the gradients the SSA/DIVA
membrane stress is assembled from (Robinson et al. 2022, Eq. 14).

Launched on `rt.grid2d`: every field involved is depth-integrated. This is the companion of
the column [`velocitygradients!`](@ref), which keeps the genuinely 3D gradients
(`z_dx`, `z_dy`, `z_dz` and the sigma-scaled `x_dz`/`y_dz`) — the two are separate kernels
because they sweep different grids, not because the work differs.

Takes no thickness argument, unlike the column method: `H` is there only for the
sigma-coordinate vertical scaling `∂/∂z = (1/H) ∂/∂ζ`, and nothing here differentiates in
z.

!!! note "Which velocity the momentum solver iterates"
    `pseudo_transient!` solves for `velocity.depthaverage_x`/`y`, not `velocity.x`/`y` —
    `ū` is genuinely 2D for SSA and DIVA alike, and leaving `velocity.x`/`y` free lets DIVA
    use them for the reconstructed 3D profile (`roadmaps/chmy.md`, Phase 3, decision 1).
"""
function depthaverage_velocitygradients!(velocity::VelocityState, rt::Runtime,
                                         mask::AbstractIceMask = NoMask())
    rt.launch2d(rt.arch, rt.grid2d,
                _depthaverage_velocity_gradients! => (velocity, mask, rt.grid2d))
    return nothing
end

# `yx`/`zx`/`zy` are the symmetric duplicates of `xy`/`xz`/`yz` and live at the same node
# classes. The collocated kernels leave them untouched; filling them costs three stores and
# removes a "why is `strainrate.yx` zero" trap for anything that reads the full tensor.
@inline function _raw_strainrate_at!(sr, vel, ::AbstractMomentumBalance, mask,
                                     I::Vararg{Integer, 3})
    if node_active(mask, NODE_AA, I[1], I[2])
        dxx = vel.x_dx[I...]
        dyy = vel.y_dy[I...]
        sr.xx[I...] = dxx
        sr.yy[I...] = dyy
        sr.zz[I...] = -(dxx + dyy)                   # incompressibility (continuity)
    else
        Z = zero(eltype(sr.xx))
        sr.xx[I...] = Z; sr.yy[I...] = Z; sr.zz[I...] = Z
    end
    _raw_strainrate_shear!(sr, vel, mask, I...)
    return nothing
end

@inline function _raw_strainrate_at!(sr, vel, ::MomentumBalance3D, mask,
                                     I::Vararg{Integer, 3})
    if node_active(mask, NODE_AA, I[1], I[2])
        sr.xx[I...] = vel.x_dx[I...]
        sr.yy[I...] = vel.y_dy[I...]
        sr.zz[I...] = vel.z_dz[I...]                 # ∂w/∂z directly
    else
        Z = zero(eltype(sr.xx))
        sr.xx[I...] = Z; sr.yy[I...] = Z; sr.zz[I...] = Z
    end
    _raw_strainrate_shear!(sr, vel, mask, I...)
    return nothing
end

@inline function _raw_strainrate_shear!(sr, vel, mask, I::Vararg{Integer, 3})
    i, j = I[1], I[2]
    Z = zero(eltype(sr.xy))
    exy = node_active(mask, NODE_AB, i, j) ?
          (vel.x_dy[I...] + vel.y_dx[I...]) / 2 : Z          # both at ab
    exz = node_active(mask, NODE_ACX_AC, i, j) ?
          (vel.x_dz[I...] + vel.z_dx[I...]) / 2 : Z          # both at acx_ac
    eyz = node_active(mask, NODE_ACY_AC, i, j) ?
          (vel.y_dz[I...] + vel.z_dy[I...]) / 2 : Z          # both at acy_ac
    sr.xy[I...] = exy; sr.yx[I...] = exy
    sr.xz[I...] = exz; sr.zx[I...] = exz
    sr.yz[I...] = eyz; sr.zy[I...] = eyz
    return nothing
end

@kernel inbounds = true function _raw_strainrate_staggered!(sr, vel, momentum, mask, O)
    I = @index(Global, NTuple)
    I = I + O
    _raw_strainrate_at!(sr, vel, momentum, mask, I...)
end

@kernel inbounds = true function _raw_strainrate_effective_staggered!(sr, mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    if node_active(mask, NODE_AA, I[1], I[2])
        # The off-diagonals are the only terms not already at `aa`, so they are the only
        # ones interpolated — once, here, rather than inside each user.
        exy = lerp(sr.xy, NODE_AA, grid, I...)
        exz = lerp(sr.xz, NODE_AA, grid, I...)
        eyz = lerp(sr.yz, NODE_AA, grid, I...)
        sr.effective[I...] = sqrt((sr.xx[I...]^2 + sr.yy[I...]^2 + sr.zz[I...]^2) / 2 +
                                  exy^2 + exz^2 + eyz^2)
    else
        sr.effective[I...] = zero(eltype(sr.effective))
    end
end

"""
$(TYPEDSIGNATURES)

Chmy-native [`raw_strainrate!`](@ref): write the strain-rate tensor components from the
velocity gradients, each at its own node class — normal components at `aa`, `ε̇_xy` at `ab`,
`ε̇_xz` at `acx`/z-`Vertex`, `ε̇_yz` at `acy`/z-`Vertex`. Every component is a sum of terms
that already live on the same node, so no interpolation enters.

`ε̇_zz` follows incompressibility (`-(ε̇_xx + ε̇_yy)`) unless `momentum` resolves the
vertical velocity ([`MomentumBalance3D`](@ref)), in which case `∂w/∂z` is used
directly — the same dispatch as the collocated method.

Does **not** compute the effective strain rate: that needs its own launch, because at `aa`
it reads `ε̇_xy` at neighbouring `ab` nodes, which the same pass may not have written yet.
Call [`raw_strainrate_effective!`](@ref) after this, or use the state-level method which
sequences all three steps.
"""
function raw_strainrate!(strainrate::StrainRateState, velocity::VelocityState,
                         momentum::AbstractMomentumBalance, rt::Runtime,
                         mask::AbstractIceMask = NoMask())
    rt.launch(rt.arch, rt.grid,
              _raw_strainrate_staggered! => (strainrate, velocity, momentum, mask))
    return nothing
end

"""
$(TYPEDSIGNATURES)

Chmy-native [`raw_strainrate_effective!`](@ref): the second invariant at `aa`, with the
off-diagonal components interpolated back from their own node classes by `lerp`.

Must run *after* [`raw_strainrate!`](@ref) — it is a stencil, not a pointwise map. Only
`interior(strainrate.effective)` is meaningful: the halo ring of the sweep reads
off-diagonal cells one further out than `raw_strainrate!`'s own sweep reached.
"""
function raw_strainrate_effective!(strainrate::StrainRateState, rt::Runtime,
                                  mask::AbstractIceMask = NoMask())
    rt.launch(rt.arch, rt.grid,
              _raw_strainrate_effective_staggered! => (strainrate, mask, rt.grid))
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level Chmy-native [`raw_strainrate!`](@ref): velocity gradients, then tensor
components, then the second invariant, in the order they depend on each other.
"""
function raw_strainrate!(mech::MechanicState, momentum::AbstractMomentumBalance,
                         rt::Runtime, mask::AbstractIceMask = NoMask())
    velocitygradients!(mech.velocity, mech.topography.thickness, rt, mask)
    raw_strainrate!(mech.strainrate, mech.velocity, momentum, rt, mask)
    raw_strainrate_effective!(mech.strainrate, rt, mask)
    return nothing
end

###############################################################
# Chmy-native, C-grid staggered membrane stress (SSA/DIVA `strainrate!`)
###############################################################
#
# This is the piece the strain-rate port above deliberately left out (see its notes in
# `roadmaps/chmy.md`, Phase 3): despite sharing the name `strainrate!` with the collocated
# dispatch it extends, the SSA/DIVA branch below is **not** the strain rate. It is the
# vertically-integrated membrane stress `2ηH·(2ε̇_xx + ε̇_yy)` and friends that the SSA/DIVA
# momentum balance actually differentiates (see the collocated method's own docstring one
# screen up). It is written into the same `StrainRateState` fields as the real strain rate
# purely for parity with the collocated code path — a legacy naming quirk kept, not fixed,
# by this port.
#
# `N_xx`, `N_yy` land at `aa`: both velocity-gradient terms they combine (`x_dx`, `y_dy`)
# already live there, so no interpolation enters, exactly like `deviatoric_stress!`'s `aa`
# branch. `N_xy` lands at `ab` and needs *two* quantities interpolated onto the corner —
# `η` harmonically (`hlerp`, the same stress-continuity argument as `deviatoric_stress!`)
# and `H` arithmetically (`lerp`, the same choice `drivingstress!` makes for thickness).
# The strict mask (`node_fully_active`) applies only to that `ab` term, for the identical
# NaN-avoidance reason as `deviatoric_stress!`: `hlerp` of `η = 0` is `NaN`, not `0`.

@kernel inbounds = true function _membrane_stress_staggered!(sxx, sxy, syy, η, H, vel,
                                                              mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    Z = zero(eltype(sxx))

    if node_active(mask, NODE_AA, i, j)
        ηH  = 2 * η[I...] * H[I...]
        dux = vel.depthaverage_x_dx[I...]
        dvy = vel.depthaverage_y_dy[I...]
        sxx[I...] = ηH * (2 * dux + dvy)
        syy[I...] = ηH * (dux + 2 * dvy)
    else
        sxx[I...] = Z
        syy[I...] = Z
    end

    sxy[I...] = node_fully_active(mask, NODE_AB, i, j) ?
        2 * hlerp(η, NODE_AB, grid, I...) * lerp(H, NODE_AB, grid, I...) *
        (vel.depthaverage_x_dy[I...] + vel.depthaverage_y_dx[I...]) / 2 : Z
end

"""
$(TYPEDSIGNATURES)

Chmy-native, C-grid staggered **membrane stress** for the SSA/DIVA momentum balance: writes
the depth-integrated components `stress.membrane_xx`/`membrane_yy` (at `aa`) and
`stress.membrane_xy` (at `ab`) from the depth-averaged viscosity
`material.viscosity_depthaveraged`, the thickness `topo.thickness` and the depth-averaged
velocity gradients already written by [`depthaverage_velocitygradients!`](@ref) — which must
run first.

The quantity is `2μ̄H·(2ε̇_xx + ε̇_yy)` and friends (Robinson et al. 2022, Eq. 14): the term
the momentum balance's stress divergence differentiates. `η` is interpolated onto `ab`
harmonically (matching [`deviatoric_stress!`](@ref) — stress, not strain rate, is continuous
across a viscosity contrast); `H` arithmetically (matching [`drivingstress!`](@ref)).

Depth-integrated throughout, so it runs on `rt.grid2d`.

!!! note "Renamed from `strainrate!`, and moved out of `StrainRateState`"
    This used to be a `strainrate!` method writing `strainrate.xx`/`xy`/`yy`, which was
    wrong twice: the quantity is a stress, not a strain rate (the collocated docstring
    admits as much), and those fields are `AA3`/`AB3` column fields, which only worked
    while `nz == 1` collapsed them onto the 2D shape this genuinely has. Both are fixed
    here (`roadmaps/chmy.md`, Phase 3, decision 2); `strainrate.xx`/`xy`/`yy` now mean only
    the true strain rate.

!!! warning "A zero viscosity gives `NaN`, not zero — same trap as `deviatoric_stress!`"
    `hlerp` averages reciprocals, so an unmasked ice-free corner produces `NaN` in
    `stress.membrane_xy` rather than `0`; pass an [`IceMask`](@ref) once the material state
    has ice-free cells.
"""
function membranestress!(stress::StressState, velocity::VelocityState,
                         material::MechanicMaterialState, topo::MechanicTopographyState,
                         momentum::MomentumBalance2D,
                         rt::Runtime, mask::AbstractIceMask = NoMask())
    rt.launch2d(rt.arch, rt.grid2d,
              _membrane_stress_staggered! =>
                  (stress.membrane_xx, stress.membrane_xy, stress.membrane_yy,
                   material.viscosity_depthaveraged, topo.thickness, velocity,
                   mask, rt.grid2d))
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level [`membranestress!`](@ref) for the SSA/DIVA momentum balance: writes
`mech.stress.membrane_xx`/`membrane_xy`/`membrane_yy` from `mech.material`,
`mech.topography` and `mech.velocity` (which must already carry the depth-averaged velocity
gradients, see [`depthaverage_velocitygradients!`](@ref)).
"""
membranestress!(mech::MechanicState,
                momentum::MomentumBalance2D,
                rt::Runtime, mask::AbstractIceMask = NoMask()) =
    membranestress!(mech.stress, mech.velocity, mech.material, mech.topography,
                    momentum, rt, mask)

###############################################################
# Chmy-native, C-grid staggered SSA/DIVA effective strain rate
###############################################################
#
# The true (not membrane-stress) second invariant ε̇_e of the *depth-averaged* velocity —
# Robinson et al. (2022), Eq. 12 — needed by `GlenViscosityContinuation`
# (`src/mechanics/solvers.jl`, wired in `src/mechanics/pseudotransient.jl`) to derive
# `viscosity_depthaveraged` from the current velocity iterate.
#
# Writes `strainrate.effective_depthaveraged` (`AA2`), *not* `strainrate.effective`
# (`AA3`). Those are two different quantities, not one at two resolutions: `.effective` is
# DIVA's Eq. 13, which carries the vertical-shear terms `¼(u_z² + v_z²)` and so genuinely
# varies layer by layer. This function computes Eq. 12, the same expression with those
# terms dropped — a single number per column. One field served both only while `nz == 1`
# collapsed `AA2` and `AA3` (`roadmaps/chmy.md`, Phase 3, decision 20).
#
# Not reusable from `raw_strainrate_effective!`: that one reads `strainrate.xy` at `ab` as
# *already written* by `raw_strainrate!`, and works on column fields throughout.

@kernel inbounds = true function _effective_strainrate_ssa_staggered!(eff, velocity, mask,
                                                                       grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    if node_active(mask, NODE_AA, i, j)
        dxx = velocity.depthaverage_x_dx[I...]
        dyy = velocity.depthaverage_y_dy[I...]
        dxy = lerp(velocity.depthaverage_x_dy, NODE_AA, grid, I...)
        dyx = lerp(velocity.depthaverage_y_dx, NODE_AA, grid, I...)
        eff[I...] = sqrt(dxx^2 + dyy^2 + dxx * dyy + ((dxy + dyx) / 2)^2)
    else
        eff[I...] = zero(eltype(eff))
    end
end

"""
$(TYPEDSIGNATURES)

Chmy-native SSA effective strain rate at `aa` (Robinson et al. 2022, Eq. 12):
`ε̇_e = √(ε̇xx² + ε̇yy² + ε̇xx·ε̇yy + ε̇xy²)`, with `ε̇xx = ∂ū/∂x`, `ε̇yy = ∂v̄/∂y` already at
`aa` and `ε̇xy = (∂ū/∂y + ∂v̄/∂x)/2` interpolated there from `ab` by `lerp`.

Writes `strainrate.effective_depthaveraged` (`AA2`) — a single value per column. DIVA's
depth-*varying* effective strain rate (Eq. 13, with the vertical-shear terms) is a
different quantity living in `strainrate.effective` (`AA3`); see the module note above.

Requires the depth-averaged velocity gradients to already be current — call
[`depthaverage_velocitygradients!`](@ref) first, exactly as [`membranestress!`](@ref) does.
"""
function effective_strainrate_ssa!(strainrate::StrainRateState, velocity::VelocityState,
                                   rt::Runtime, mask::AbstractIceMask = NoMask())
    rt.launch2d(rt.arch, rt.grid2d,
              _effective_strainrate_ssa_staggered! =>
                  (strainrate.effective_depthaveraged, velocity, mask, rt.grid2d))
    return nothing
end
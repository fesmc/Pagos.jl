"""
$(TYPEDSIGNATURES)

Compute the strain rate tensor in place.
"""
function strainrate!(
    strainrate,
    velocity,
    material,
    topo,
    momentum::AbstractMomentumBalance,
)
    backend = get_backend(strainrate.xx)
    kernel! = _strainrate_kernel!(backend)
    kernel!(
        strainrate,
        velocity,
        material,
        topo,
        momentum;
        ndrange = length(strainrate.xx),
    )
    return nothing
end

# Per-element fallback: triggers when a momentum balance has no specialized method below,
# i.e. the real extension point. Guarding the launcher instead would never catch this.
function strainrate!(
    strainrate,
    velocity,
    material,
    topo,
    momentum::AbstractMomentumBalance,
    I,
)
    throw(ArgumentError("Unsupported momentum balance type: $(typeof(momentum))"))
end

function strainrate!(
    strainrate,
    velocity,
    material,
    topo,
    momentum::SIAMomentumBalance,
    I,
)
    # For SIA: depth-averaged viscosity, thickness, strainrate.xz/yz, velocity.*_bar_dz are all 2D
    strainrate.xz[I] =
        material.viscosity_depthaveraged[I] *
        topo.thickness[I] *
        velocity.depthaverage_x_dz[I]
    strainrate.yz[I] =
        material.viscosity_depthaveraged[I] *
        topo.thickness[I] *
        velocity.depthaverage_y_dz[I]
end

function strainrate!(
    strainrate,
    velocity,
    material,
    topo,
    momentum::MomentumBalance2D,
    I,
)
    # For SSA/DIVA: depth-averaged viscosity, thickness, strainrate.xx/xy/yy, velocity.*_d* are all 2D
    strainrate.xx[I] =
        2 *
        material.viscosity_depthaveraged[I] *
        topo.thickness[I] *
        (2 * velocity.x_dx[I] + velocity.y_dy[I])
    strainrate.xy[I] =
        material.viscosity_depthaveraged[I] *
        topo.thickness[I] *
        (velocity.x_dy[I] + velocity.y_dx[I])
    strainrate.yy[I] =
        2 *
        material.viscosity_depthaveraged[I] *
        topo.thickness[I] *
        (velocity.x_dx[I] + 2 * velocity.y_dy[I])
end

function strainrate!(
    strainrate,
    velocity,
    material,
    topo,
    momentum::BlatterPattynMomentumBalance,
    I,
)
    # For Blatter-Pattyn: 3D viscosity, strainrate.xx/xy/yy/xz/yz, velocity.*_d* are all 3D
    strainrate.xx[I] =
        2 * material.viscosity[I] * (2 * velocity.x_dx[I] + velocity.y_dy[I])
    strainrate.xy[I] = material.viscosity[I] * (velocity.x_dy[I] + velocity.y_dx[I])
    strainrate.yy[I] =
        2 * material.viscosity[I] * (velocity.x_dx[I] + 2 * velocity.y_dy[I])
    strainrate.xz[I] = material.viscosity[I] * velocity.x_dz[I]
    strainrate.yz[I] = material.viscosity[I] * velocity.y_dz[I]
end

@kernel function _strainrate_kernel!(
    strainrate,
    velocity,
    material,
    topo,
    momentum::AbstractMomentumBalance,
)
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
function strainrate_effective!(
    strainrate,
    velocity,
    momentum::AbstractMomentumBalance,
    I,
)
    throw(ArgumentError("Unsupported momentum balance type: $(typeof(momentum))"))
end

function strainrate_effective!(strainrate, velocity, momentum::SIAMomentumBalance, I)
    strainrate.effective[I] = sqrt(
        1 / 4 * (velocity.depthaverage_x_dz[I] + velocity.depthaverage_y_dz[I]) ^ 2,
    )
end

function strainrate_effective!(strainrate, velocity, momentum::SSAMomentumBalance, I)
    strainrate.effective[I] = sqrt(
        velocity.x_dx[I]^2 +
        velocity.y_dy[I]^2 +
        velocity.x_dx[I] * velocity.y_dy[I] +
        1 / 4 * (velocity.x_dy[I] + velocity.y_dx[I]) ^ 2,
    )
end

# Deliberately a custom union, not [`MomentumBalance2D`](@ref)/[`MomentumBalance3D`](@ref):
# groups by whether ε̇_e includes the vertical-shear terms ¼(u_z² + v_z²) (DIVA and BP share
# this, Robinson et al. 2022 Eq. 13) rather than by unknown dimensionality (SSA doesn't,
# despite being in the same 2D group as DIVA).
function strainrate_effective!(
    strainrate,
    velocity,
    momentum::MB,
    I,
) where {MB<:Union{DIVAMomentumBalance,BlatterPattynMomentumBalance}}
    strainrate.effective[I] = sqrt(
        velocity.x_dx[I]^2 +
        velocity.y_dy[I]^2 +
        velocity.x_dx[I] * velocity.y_dy[I] +
        1 / 4 * (velocity.x_dy[I] + velocity.y_dx[I]) ^ 2 +
        1 / 4 * velocity.x_dz[I]^2 +
        1 / 4 * velocity.y_dz[I]^2,
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

@kernel function _raw_strainrate_kernel!(
    strainrate,
    velocity,
    momentum::AbstractMomentumBalance,
)
    I = @index(Global, Linear)
    @inbounds begin
        raw_strainrate!(strainrate, velocity, momentum, I)
        raw_strainrate_effective!(strainrate, momentum, I)
    end
end

# Default (plane / depth-integrated): ε̇_zz reconstructed from incompressibility.
#
# `I::Integer` (the flat `@index(Global, Linear)` from the kernel above) is needed only to
# disambiguate this 4-argument per-element method from the 4-argument staggered
# `raw_strainrate!(sr, vel, momentum, rt::Runtime)` further down — left untyped, the two
# would match the same call with neither more specific.
function raw_strainrate!(
    strainrate,
    velocity,
    momentum::AbstractMomentumBalance,
    I::Integer,
)
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
function raw_strainrate!(strainrate, velocity, momentum::MomentumBalance3D, I::Integer)
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
        strainrate.xy[I]^2 +
        strainrate.xz[I]^2 +
        strainrate.yz[I]^2,
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
    return velocitygradients!(
        velocity.x_dx,
        velocity.x_dy,
        velocity.y_dx,
        velocity.y_dy,
        x,
        y,
        dx,
        dy,
    )
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
# The point of the C-grid layout: assigning `u` to `acx`, `v` to `acy` and `w` to
# `aa`/z-`Vertex` fixes every gradient's node by operator algebra, so **no interpolation
# appears anywhere in the tensor**:
#
#   ∂u/∂x : acx    → aa       ∂v/∂x : acy    → ab       ∂w/∂x : aa_ac → acx_ac
#   ∂u/∂y : acx    → ab       ∂v/∂y : acy    → aa       ∂w/∂y : aa_ac → acy_ac
#   ∂u/∂z : acx    → acx_ac   ∂v/∂z : acy    → acy_ac   ∂w/∂z : aa_ac → aa
#
# Each strain-rate component is then a sum of terms already on the *same* node (e.g.
# ε̇_xy = (∂u/∂y + ∂v/∂x)/2, both at `ab`). The interpolations that remain are exactly where
# physics genuinely mixes node classes: viscosity (lives at `aa`, scales off-diagonal
# strain rates elsewhere) and the second invariants (cell-centred, built from off-diagonal
# components that are not).
#
# Two things this forces, versus the collocated kernels above:
#
#  1. No shared flat `@index(Global, Linear)` — the four node classes (`aa`, `ab`, `acx_ac`,
#     `acy_ac`) have different shapes, so each kernel indexes every field at its own
#     `(i, j, k)` instead. Still one fused kernel; only the flat index goes.
#  2. The invariants need their own launch: `ε̇_e` at `aa` reads `ε̇_xy` at neighbouring `ab`
#     nodes, which may not be written yet in the same pass. Hence `raw_strainrate!` then
#     `raw_strainrate_effective!`, in that order.

# Physical vertical derivative on the terrain-following sigma axis: `∂/∂z = (1/H) ∂/∂ζ`.
# `∂z_σ` supplies `∂/∂ζ` — Chmy's own `∂z` is wrong on a non-uniform axis, see
# `src/api/sigma_operators.jl`. Returns zero rather than `Inf`/`NaN` where there is no ice,
# since an `Inf` here would propagate into the whole tensor.
@inline _dz_over_H(dζ, H) = H > zero(H) ? dζ / H : zero(dζ)

@kernel inbounds = true function _velocity_gradients!(
    velocity,
    H,
    mask,
    grid,
    grid2d,
    O,
)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    u, v, w = velocity.x, velocity.y, velocity.z
    Z = zero(eltype(velocity.x_dx))

    # `node_active` only consults the horizontal part of the node class, so the `_AC`
    # (z-`Vertex`) variants share activity with `NODE_ACX`/`NODE_ACY`.
    act_aa = node_active(mask, NODE_AA, i, j)
    act_ab = node_active(mask, NODE_AB, i, j)
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
    velocity.x_dz[I...] =
        act_acx ? _dz_over_H(∂z_σ(u, grid, I...), lerp(H, NODE_ACX, grid2d, i, j, 1)) :
        Z
    velocity.y_dz[I...] =
        act_acy ? _dz_over_H(∂z_σ(v, grid, I...), lerp(H, NODE_ACY, grid2d, i, j, 1)) :
        Z
    velocity.z_dz[I...] = act_aa ? _dz_over_H(∂z_σ(w, grid, I...), H[i, j, 1]) : Z
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

!!! warning "Horizontal derivatives are at constant ζ until [`terrain_metric_correction!`](@ref) runs"
    `∂u/∂x` here is taken at constant ζ, not at constant z. On the
    [`MomentumBalance3D`](@ref) path [`pseudo_rate!`](@ref) applies the terrain-following
    correction immediately afterwards, in a second pass (it cannot be fused — see that
    function's docstring). SSA/DIVA never call it: their gradients are depth-averaged and
    carry no `ζ` dependence to correct.
"""
function velocitygradients!(
    velocity::VelocityState,
    H,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    rt.launch(
        rt.arch,
        rt.grid,
        _velocity_gradients! => (velocity, H, mask, rt.grid, rt.grid2d),
    )
    return nothing
end

###############################################################
# Terrain-following (sigma) metric correction on the horizontal gradients
###############################################################
#
# `_velocity_gradients!` above differences at constant ζ; every balance wants the derivative
# at constant *z*. For the terrain-following map `z = b + ζH = s - (1-ζ)H` the two differ by
# one term (and the mirror image in y):
#
#   ∂f/∂x|_z = ∂f/∂x|_ζ - c_x ∂f/∂z,     c_x ≡ ∂z/∂x|_ζ = ∂s/∂x - (1-ζ) ∂H/∂x
#
# `c_x` is the local slope of the ζ-surface the difference was taken along, zero only where
# the surface and thickness are both flat. Not a small correction on the geometries BP
# exists for: on ISMIP-HOM B at `L = 10 km`, `∂H/∂x` reaches 0.31, so `c_x ∂u/∂z` exceeds the
# uncorrected `∂u/∂x|_ζ`. Dropping it is defensible for SSA/DIVA (shallow by construction),
# not for BP.
#
# Cannot be fused into `_velocity_gradients!`: the correction to `∂u/∂x` at `aa` reads
# `∂u/∂z` at the four surrounding `acx_ac` nodes, which that kernel is still writing in the
# same sweep. Hence a second pass, same argument as `raw_strainrate_effective!`.
#
# `c_x` is evaluated from `s` and `H` (not a bed field) because those are what
# `MechanicTopographyState` carries, and `b = s - H` makes the identity exact.
#
# @dev: only the *inner* gradient is corrected here. The **outer** divergence
# `∂σxx/∂x|_z = ∂x(σxx)|_ζ - c_x ∂σxx/∂z` carries the same term and is not corrected yet
# (`_dotvel_staggered_bp!` still differences at constant ζ). Conservative form:
#
#   r_x = (1/H)[ ∂(H σxx)/∂x|_ζ + ∂(H σxy)/∂y|_ζ + ∂(σxz - σxx c_x - σxy c_y)/∂ζ ]
#
# — the outer correction folds into the vertical flux, which also makes the stress-free
# surface condition the true `τ·n = 0` (equations doc, item A2). That form needs `σxx`/`σxy`
# at the layer interfaces and a matching Gershgorin row sum, so it's left for its own change:
# an unbounded new coupling in the operator is more than the explicit iteration can absorb
# right now.

@kernel inbounds = true function _terrain_metric_correction!(
    velocity,
    H,
    s,
    mask,
    grid,
    grid2d,
    O,
)
    I = @index(Global, NTuple)
    I = I + O
    i, j, k = I
    T = eltype(velocity.x_dx)
    w = one(T) - T(zcenter(grid, k))          # (1 - ζ) at the layer midpoint

    ## Written out rather than delegated to `lerp`: `x_dz`/`y_dz` are staggered in *two* axes
    ## at once relative to their targets, and Chmy's `itp` does not unroll that case under
    ## `GPUCompiler` (dynamic tuple index inside `ntuple` → `InvalidIRError` on a `CuArray`).
    ## Plain 4-point average, pinned by the uniform-slab and staggered-gradient tests.
    q4(f, a, b, c, d) = (f[a...] + f[b...] + f[c...] + f[d...]) / 4

    ## This kernel reads one node further out than `velocitygradients!` writes, so on the
    ## outermost launched ring it folds in entries nothing ever wrote (a deterministic zero),
    ## halving the correction there; contamination stays in the halo, not the interior.
    ## @dev: do not "fix" this by clamping or skipping — both were tried and both are worse
    ## (clamping imports `NaN` and broke the masked `ImplicitVertical` fixed-point test;
    ## skipping the ring broke the uniform-slab test). The real fix is a proper halo/BC
    ## treatment for the gradient fields (`pagos-roadmap/blatter-pattyn-equations.md`, A1 and A6).

    if node_active(mask, NODE_AA, i, j)
        ## `∂x` of an `aa` field lands at `acx`; averaging the two faces of the cell is the
        ## centred difference at `aa`.
        cx =
            (∂x(s, grid2d, i, j, 1) + ∂x(s, grid2d, i + 1, j, 1)) / 2 -
            w * (∂x(H, grid2d, i, j, 1) + ∂x(H, grid2d, i + 1, j, 1)) / 2
        cy =
            (∂y(s, grid2d, i, j, 1) + ∂y(s, grid2d, i, j + 1, 1)) / 2 -
            w * (∂y(H, grid2d, i, j, 1) + ∂y(H, grid2d, i, j + 1, 1)) / 2
        ## acx_ac → aa: average the two `x` faces and the two `ζ` interfaces of the cell.
        uz = q4(
            velocity.x_dz,
            (i, j, k),
            (i + 1, j, k),
            (i, j, k + 1),
            (i + 1, j, k + 1),
        )
        ## acy_ac → aa: the same in `y`.
        vz = q4(
            velocity.y_dz,
            (i, j, k),
            (i, j + 1, k),
            (i, j, k + 1),
            (i, j + 1, k + 1),
        )
        velocity.x_dx[I...] -= cx * uz
        velocity.y_dy[I...] -= cy * vz
    end

    if node_active(mask, NODE_AB, i, j)
        ## At `ab` the same two gradients are needed one node over: `∂x(s)` is at `acx`, so
        ## it is averaged in *y*; `∂y(s)` is at `acy`, so it is averaged in *x*.
        cx =
            (∂x(s, grid2d, i, j - 1, 1) + ∂x(s, grid2d, i, j, 1)) / 2 -
            w * (∂x(H, grid2d, i, j - 1, 1) + ∂x(H, grid2d, i, j, 1)) / 2
        cy =
            (∂y(s, grid2d, i - 1, j, 1) + ∂y(s, grid2d, i, j, 1)) / 2 -
            w * (∂y(H, grid2d, i - 1, j, 1) + ∂y(H, grid2d, i, j, 1)) / 2
        ## acx_ac → ab: `x` already sits on the vertex, so only `y` and `ζ` are averaged.
        uz = q4(
            velocity.x_dz,
            (i, j - 1, k),
            (i, j, k),
            (i, j - 1, k + 1),
            (i, j, k + 1),
        )
        ## acy_ac → ab: `y` already on the vertex; average `x` and `ζ`.
        vz = q4(
            velocity.y_dz,
            (i - 1, j, k),
            (i, j, k),
            (i - 1, j, k + 1),
            (i, j, k + 1),
        )
        velocity.x_dy[I...] -= cy * uz
        velocity.y_dx[I...] -= cx * vz
    end
end

"""
$(TYPEDSIGNATURES)

Convert the four horizontal velocity gradients from constant-ζ to constant-`z` derivatives,
in place — the terrain-following metric correction
`∂f/∂x|_z = ∂f/∂x|_ζ - (∂s/∂x - (1-ζ)∂H/∂x)·∂f/∂z` and its `y` mirror
(`pagos-roadmap/blatter-pattyn-equations.md`, item A1).

Corrects `x_dx`/`y_dy` (at `aa`) and `x_dy`/`y_dx` (at `ab`); `x_dz`/`y_dz` are already
physical `∂/∂z` and are read, not written. Must run **after** [`velocitygradients!`](@ref)
and before anything reading the gradients — it is a separate launch because it reads
`x_dz`/`y_dz` at neighbouring nodes that the gradient kernel is still writing.

Identically zero wherever the surface and thickness are both laterally uniform, so it does
not move a uniform-slab result.

!!! note "The inner gradients only"
    The matching correction to the *outer* stress divergence is not applied — see the source
    note above for the conservative form it takes and why it is a separate change.
"""
function terrain_metric_correction!(
    velocity::VelocityState,
    H,
    s,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    rt.launch(
        rt.arch,
        rt.grid,
        _terrain_metric_correction! => (velocity, H, s, mask, rt.grid, rt.grid2d),
    )
    return nothing
end

###############################################################
# Velocity-gradient clamp — a Blatter-Pattyn numerical safety net
###############################################################
#
# Not a physical strain-rate regularization (contrast `GlenViscosityContinuation`'s `ε̇0`,
# which floors a *denominator*): a hard ceiling on the gradients themselves, for real
# geometry where a masked, ice-free-adjacent column can lose its entire vertical-stiffness
# contribution to the Gershgorin bound (`node_fully_active` zeroes `R₋`/`R₊` at every `k` in
# that column, found on real 8 km AIS geometry). The vertical term is normally dominant, so
# losing it can leave `Δτ` orders of magnitude too large there, and a single explicit step
# produces an unphysical velocity — which the next iteration reads back as an unphysical
# gradient, feeding an equally unphysical stress right back into the residual.
#
# Clamping the gradients breaks that feedback loop at its source, complementing (not
# replacing) `gershgorin_dt!`'s own `dtau_cap`, which addresses the first step's `Δτ`
# directly. Real ice strain rates are `~1e-4`–`~1e-2` /yr even in fast shear margins, so any
# `cap` worth using sits far below where this could ever bind on a physically sane velocity
# field; `cap = Inf` (the default everywhere this is threaded through) is a no-op.

@kernel inbounds = true function _clamp_velocity_gradients!(velocity, cap, O)
    I = @index(Global, NTuple)
    I = I + O
    velocity.x_dx[I...] = clamp(velocity.x_dx[I...], -cap, cap)
    velocity.y_dy[I...] = clamp(velocity.y_dy[I...], -cap, cap)
    velocity.x_dy[I...] = clamp(velocity.x_dy[I...], -cap, cap)
    velocity.y_dx[I...] = clamp(velocity.y_dx[I...], -cap, cap)
    velocity.x_dz[I...] = clamp(velocity.x_dz[I...], -cap, cap)
    velocity.y_dz[I...] = clamp(velocity.y_dz[I...], -cap, cap)
end

"""
$(TYPEDSIGNATURES)

Clamp every velocity-gradient component (`x_dx`, `y_dy`, `x_dy`, `y_dx`, `x_dz`, `y_dz`) to
`[-cap, cap]`, in place. A [`MomentumBalance3D`](@ref) numerical safety net, not a physical
regularization — see the source note above. `cap = Inf` (the default) is a no-op; masked
(inactive) nodes are already zero from [`velocitygradients!`](@ref), so clamping them is
harmless and the mask is not re-checked here.

Call after [`velocitygradients!`](@ref) and before anything that reads the gradients
(membrane stress, the effective strain rate) — [`pseudo_rate!`](@ref) does so for BP.
"""
function clamp_velocity_gradients!(velocity::VelocityState, cap, rt::Runtime)
    isfinite(cap) || return nothing
    rt.launch(rt.arch, rt.grid, _clamp_velocity_gradients! => (velocity, cap))
    return nothing
end

# The depth-averaged counterpart of `_velocity_gradients!`: same node algebra (`∂x` of an
# `acx` field lands on `aa`, `∂y` on `ab`), but every field is depth-integrated, so this runs
# on `grid2d` with no thickness argument — the sigma scaling `H` exists for in the column
# kernel has no counterpart here.
@kernel inbounds = true function _depthaverage_velocity_gradients!(
    velocity,
    mask,
    grid,
    O,
)
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
    use them for the reconstructed 3D profile.
"""
function depthaverage_velocitygradients!(
    velocity::VelocityState,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _depthaverage_velocity_gradients! => (velocity, mask, rt.grid2d),
    )
    return nothing
end

###############################################################
# Blatter-Pattyn/Stokes vertical velocity from incompressibility
###############################################################
#
# `w` is not an unknown of the momentum solve: BP resolves `u`/`v` only, and
# `∂w/∂z = -(u_x + v_y)` diagnoses `w` afterward by integrating up from the bed where
# `w = 0` (no basal melting/penetration). Per-column, serial in `k` — same `rt.launch2d` +
# internal `for k in 1:nz` shape as `_viscosity_integrals!` (`src/mechanics/velocities.jl`),
# for the same reason: a cumulative sum has no useful parallelism across `k`.
#
# Midpoint rule: `Δw_k = -(u_x[k] + v_y[k]) · Δz_k(phys)`, `Δz_k(phys) = Δζ_k · H` — the same
# sigma-to-physical conversion `_viscosity_integrals!` uses for its own `dz`. Requires
# `velocity.x_dx`/`y_dy` to already be current — call [`velocitygradients!`](@ref) first.

@kernel inbounds = true function _verticalvelocity!(w, dux, dvy, H, nz, mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    T = eltype(w)
    Hij = H[i, j, 1]

    if node_active(mask, NODE_AA, i, j) && Hij > zero(T)
        acc = zero(T)
        w[i, j, 1] = acc
        for k = 1:nz
            acc -= (dux[i, j, k] + dvy[i, j, k]) * Δz(grid, Center(), i, j, k) * Hij
            w[i, j, k+1] = acc
        end
    else
        for k = 1:(nz+1)
            w[i, j, k] = zero(T)
        end
    end
end

"""
$(TYPEDSIGNATURES)

Diagnose the vertical velocity `velocity.z` (`AAZ3`, at the layer interfaces) from
incompressibility, `∂w/∂z = -(u_x + v_y)`, integrated up from the bed (`w = 0`) —
[`MomentumBalance3D`](@ref)'s `w` is a post-solve diagnostic, not an unknown of
[`pseudo_transient!`](@ref) (`pagos-roadmap/blatter-pattyn.md`, Phase 1), so this is called
*after* a converged solve, never from inside the PT loop.

Requires `velocity.x_dx`/`y_dy` to already be current (see [`velocitygradients!`](@ref)).
Feeds `velocity.z_dx`/`z_dy`/`z_dz` on the *next* call to [`velocitygradients!`](@ref) — that
function already differentiates `velocity.z` (see `_velocity_gradients!`'s `z_dx`/`z_dy`/
`z_dz` lines), so nothing here duplicates that; it exists only to fill the `w` those read.
"""
function verticalvelocity!(
    velocity::VelocityState,
    H,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    nz = size(rt.grid, Center())[3]
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _verticalvelocity! =>
            (velocity.z, velocity.x_dx, velocity.y_dy, H, nz, mask, rt.grid),
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level [`verticalvelocity!`](@ref).
"""
verticalvelocity!(mech::MechanicState, rt::Runtime, mask::AbstractIceMask = NoMask()) =
    verticalvelocity!(mech.velocity, mech.topography.thickness, rt, mask)

@inline function _raw_strainrate_at!(
    sr,
    vel,
    ::AbstractMomentumBalance,
    mask,
    I::Vararg{Integer,3},
)
    if node_active(mask, NODE_AA, I[1], I[2])
        dxx = vel.x_dx[I...]
        dyy = vel.y_dy[I...]
        sr.xx[I...] = dxx
        sr.yy[I...] = dyy
        sr.zz[I...] = -(dxx + dyy)                   # incompressibility (continuity)
    else
        Z = zero(eltype(sr.xx))
        sr.xx[I...] = Z
        sr.yy[I...] = Z
        sr.zz[I...] = Z
    end
    _raw_strainrate_shear!(sr, vel, mask, I...)
    return nothing
end

@inline function _raw_strainrate_at!(
    sr,
    vel,
    ::MomentumBalance3D,
    mask,
    I::Vararg{Integer,3},
)
    if node_active(mask, NODE_AA, I[1], I[2])
        sr.xx[I...] = vel.x_dx[I...]
        sr.yy[I...] = vel.y_dy[I...]
        sr.zz[I...] = vel.z_dz[I...]                 # ∂w/∂z directly
    else
        Z = zero(eltype(sr.xx))
        sr.xx[I...] = Z
        sr.yy[I...] = Z
        sr.zz[I...] = Z
    end
    _raw_strainrate_shear!(sr, vel, mask, I...)
    return nothing
end

@inline function _raw_strainrate_shear!(sr, vel, mask, I::Vararg{Integer,3})
    i, j = I[1], I[2]
    Z = zero(eltype(sr.xy))
    exy = node_active(mask, NODE_AB, i, j) ? (vel.x_dy[I...] + vel.y_dx[I...]) / 2 : Z          # both at ab
    exz =
        node_active(mask, NODE_ACX_AC, i, j) ? (vel.x_dz[I...] + vel.z_dx[I...]) / 2 : Z          # both at acx_ac
    eyz =
        node_active(mask, NODE_ACY_AC, i, j) ? (vel.y_dz[I...] + vel.z_dy[I...]) / 2 : Z          # both at acy_ac
    sr.xy[I...] = exy
    sr.xz[I...] = exz
    sr.yz[I...] = eyz
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
        sr.effective[I...] = sqrt(
            (sr.xx[I...]^2 + sr.yy[I...]^2 + sr.zz[I...]^2) / 2 + exy^2 + exz^2 + eyz^2,
        )
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
function raw_strainrate!(
    strainrate::StrainRateState,
    velocity::VelocityState,
    momentum::AbstractMomentumBalance,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    rt.launch(
        rt.arch,
        rt.grid,
        _raw_strainrate_staggered! => (strainrate, velocity, momentum, mask),
    )
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
function raw_strainrate_effective!(
    strainrate::StrainRateState,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    rt.launch(
        rt.arch,
        rt.grid,
        _raw_strainrate_effective_staggered! => (strainrate, mask, rt.grid),
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level Chmy-native [`raw_strainrate!`](@ref): velocity gradients, then tensor
components, then the second invariant, in the order they depend on each other.
"""
function raw_strainrate!(
    mech::MechanicState,
    momentum::AbstractMomentumBalance,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    velocitygradients!(mech.velocity, mech.topography.thickness, rt, mask)
    raw_strainrate!(mech.strainrate, mech.velocity, momentum, rt, mask)
    raw_strainrate_effective!(mech.strainrate, rt, mask)
    return nothing
end

###############################################################
# Chmy-native, C-grid staggered membrane stress (SSA/DIVA)
###############################################################
#
# Despite sharing the name `membranestress!` with the collocated dispatch it extends, this
# is **not** the strain rate: it's the vertically-integrated membrane stress
# `2ηH·(2ε̇_xx + ε̇_yy)` and friends that the SSA/DIVA momentum balance actually
# differentiates (see the collocated method's own docstring one screen up).
#
# `N_xx`, `N_yy` land at `aa`: both velocity-gradient terms they combine (`x_dx`, `y_dy`)
# already live there, so no interpolation enters, exactly like `deviatoric_stress!`'s `aa`
# branch. `N_xy` lands at `ab` and needs *two* quantities interpolated onto the corner —
# `η` harmonically (`hlerp`, the same stress-continuity argument as `deviatoric_stress!`)
# and `H` arithmetically (`lerp`, the same choice `drivingstress!` makes for thickness).
# The strict mask (`node_fully_active`) applies only to that `ab` term, for the identical
# NaN-avoidance reason as `deviatoric_stress!`: `hlerp` of `η = 0` is `NaN`, not `0`.
#
# The three `@inline` helpers below exist so the *same* algebra serves both the
# self-contained kernel (`_membrane_stress_staggered!`, which recomputes the `η`/`H`
# prefactor every call) and the cached pair (`_membrane_prefactors!` +
# `_membrane_stress_cached!`, which the PT loop uses — see `membrane_prefactors!`). Writing
# the split form out a second time by hand is exactly how a factor of two goes missing: note
# that the `ab` prefactor is `hlerp·lerp`, *not* `2·hlerp·lerp`, because the `2` and the `/2`
# in `N_xy = 2·η̄·H̄·(u_y + v_x)/2` cancel. Multiplying and dividing by 2 are exact in binary
# floating point, so the cached path is bit-for-bit identical to the direct one.

# `2ηH` at `aa` — the prefactor of both normal components.
@inline _membrane_pre_aa(η, H, mask, I::Vararg{Integer,N}) where {N} =
    node_active(mask, NODE_AA, I[1], I[2]) ? 2 * η[I...] * H[I...] : zero(eltype(η))

# `η̄H̄` at `ab`, `η` harmonically and `H` arithmetically interpolated onto the corner.
@inline _membrane_pre_ab(η, H, mask, grid, I::Vararg{Integer,N}) where {N} =
    node_fully_active(mask, NODE_AB, I[1], I[2]) ?
    hlerp(η, NODE_AB, grid, I...) * lerp(H, NODE_AB, grid, I...) : zero(eltype(η))

# The velocity-gradient half, given both prefactors already evaluated at `I`.
#
# The masks are re-tested here rather than left to `p_aa`/`p_ab` being zero, and that is not
# redundant: `0.0 * g` is `-0.0` for any `g < 0`, while the branch writes `+0.0`. Both are
# arithmetically zero and nothing downstream can tell them apart, but reproducing the
# original's exact bit pattern is what lets this path be swapped in under a bit-for-bit test
# rather than an approximate one. The mask reads are Bool loads already in cache from the
# rest of the loop; the cost this hoists out is `hlerp`'s divisions, not these.
@inline function _membrane_from_pre(p_aa, p_ab, vel, mask, I::Vararg{Integer,N}) where {N}
    Z = zero(p_aa)
    if node_active(mask, NODE_AA, I[1], I[2])
        dux = vel.depthaverage_x_dx[I...]
        dvy = vel.depthaverage_y_dy[I...]
        sxx = p_aa * (2 * dux + dvy)
        syy = p_aa * (dux + 2 * dvy)
    else
        sxx = Z
        syy = Z
    end
    sxy =
        node_fully_active(mask, NODE_AB, I[1], I[2]) ?
        p_ab * (vel.depthaverage_x_dy[I...] + vel.depthaverage_y_dx[I...]) : Z
    return sxx, sxy, syy
end

@kernel inbounds = true function _membrane_stress_staggered!(
    sxx,
    sxy,
    syy,
    η,
    H,
    vel,
    mask,
    grid,
    O,
)
    I = @index(Global, NTuple)
    I = I + O
    a, b, c = _membrane_from_pre(
        _membrane_pre_aa(η, H, mask, I...),
        _membrane_pre_ab(η, H, mask, grid, I...),
        vel,
        mask,
        I...,
    )
    sxx[I...] = a
    sxy[I...] = b
    syy[I...] = c
end

# Writes the two prefactors the PT loop caches. Same expressions as the kernel above reads
# inline, via the same helpers — see `membrane_prefactors!` for when this is refreshed.
@kernel inbounds = true function _membrane_prefactors!(p_aa, p_ab, η, H, mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    p_aa[I...] = _membrane_pre_aa(η, H, mask, I...)
    p_ab[I...] = _membrane_pre_ab(η, H, mask, grid, I...)
end

# The hot-loop kernel: no `hlerp` (nine FP64 divisions a node), no `lerp`, no `η`/`H` reads
# — just two field loads and the velocity-gradient arithmetic.
@kernel inbounds = true function _membrane_stress_cached!(
    sxx,
    sxy,
    syy,
    p_aa,
    p_ab,
    vel,
    mask,
    O,
)
    I = @index(Global, NTuple)
    I = I + O
    a, b, c = _membrane_from_pre(p_aa[I...], p_ab[I...], vel, mask, I...)
    sxx[I...] = a
    sxy[I...] = b
    syy[I...] = c
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

!!! warning "A zero viscosity gives `NaN`, not zero — same trap as `deviatoric_stress!`"
    `hlerp` averages reciprocals, so an unmasked ice-free corner produces `NaN` in
    `stress.membrane_xy` rather than `0`; pass an [`IceMask`](@ref) once the material state
    has ice-free cells.
"""
function membranestress!(
    stress::StressState,
    velocity::VelocityState,
    material::MechanicMaterialState,
    topo::MechanicTopographyState,
    momentum::MomentumBalance2D,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _membrane_stress_staggered! => (
            stress.membrane_xx,
            stress.membrane_xy,
            stress.membrane_yy,
            material.viscosity_depthaveraged,
            topo.thickness,
            velocity,
            mask,
            rt.grid2d,
        ),
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level [`membranestress!`](@ref) for the SSA/DIVA momentum balance: writes
`mech.stress.membrane_xx`/`membrane_xy`/`membrane_yy` from `mech.material`,
`mech.topography` and `mech.velocity` (which must already carry the depth-averaged velocity
gradients, see [`depthaverage_velocitygradients!`](@ref)).
"""
membranestress!(
    mech::MechanicState,
    momentum::MomentumBalance2D,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
) = membranestress!(
    mech.stress,
    mech.velocity,
    mech.material,
    mech.topography,
    momentum,
    rt,
    mask,
)

"""
$(TYPEDSIGNATURES)

Write the cached membrane prefactors `2η̄H` (at `aa`) and `η̄H̄` (at `ab`) into
`solver.membrane_pre_aa`/`membrane_pre_ab`, from `material.viscosity_depthaveraged` and
`topo.thickness`.

This is the expensive half of [`membranestress!`](@ref) — `hlerp` alone is nine FP64
divisions per node, and at Float64 the membrane kernel is the PT loop's single most costly
one (32% of GPU device time, `benchmark/basics/gpu/README.md` §3). Neither `µ̄` nor `H` moves
during most solves, so the loop builds this once and then runs the cheap cached
[`membranestress!`](@ref) method every iteration.

Measured on the whole `pseudo_transient!` loop (CPU, 4 threads, 380×380×4, Float64, 100
fixed iterations, best of 5), against recomputing the prefactor every iteration:

| configuration | before | after | |
|---|---:|---:|---:|
| DIVA + `NoDIVUpdate` (`µ̄` static) | 10.93 | 7.84 ms/iter | **1.39×** |
| SSA + `NoViscosityContinuation` (`µ̄` static) | 9.52 | 7.79 ms/iter | **1.22×** |
| DIVA + `PeriodicDIVUpdate(10)` | 15.29 | 13.11 ms/iter | 1.17× |
| SSA + `GlenViscosityContinuation` (`µ̄` every iteration) | 15.01 | 13.81 ms/iter | 1.09× |

The last two rebuild the cache as often as `µ̄` moves and so cannot benefit from the hoist
itself; they still come out ahead because splitting the work leaves two simpler kernels than
the one fused kernel they replace. There is therefore no configuration this costs, and no
need for a second dispatch path that skips the cache.

!!! warning "Whoever writes `µ̄` owns refreshing this"
    The cache is only as good as its invalidation. `H` never moves inside a momentum solve,
    but `µ̄` is written by [`update_viscosity!`](@ref) — every iteration under
    [`GlenViscosityContinuation`](@ref), and at [`PeriodicDIVUpdate`](@ref)'s cadence inside
    [`diva_update!`](@ref). [`pseudo_transient!`](@ref) calls
    [`_refresh_membrane_prefactors!`](@ref) after both, which dispatches to a no-op when the
    continuation is [`NoViscosityContinuation`](@ref) and so cannot have moved `µ̄`. A new
    writer of `µ̄` must do the same or the solve silently uses a stale viscosity.
"""
function membrane_prefactors!(
    solver::PseudoTransientSolver,
    mech::MechanicState,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _membrane_prefactors! => (
            solver.membrane_pre_aa,
            solver.membrane_pre_ab,
            mech.material.viscosity_depthaveraged,
            mech.topography.thickness,
            mask,
            rt.grid2d,
        ),
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

[`membranestress!`](@ref) from the prefactors [`membrane_prefactors!`](@ref) already wrote,
rather than recomputing `2η̄H`/`η̄H̄` from `µ̄` and `H`.

**Bit-for-bit identical** to the self-contained method, and verified as such over a full
`pseudo_transient!` solve in every viscosity-continuation regime — the two share
`_membrane_from_pre`, the `2`/`/2` that move between the halves are exact in binary floating
point, and the mask is re-tested here so masked nodes get `+0.0` rather than the `-0.0` that
`0.0 * negative` would produce.

`prefactors` is the `solver`, passed for its two cache fields; the argument exists so this
method cannot be reached by a caller that has not built them.
"""
function membranestress!(
    stress::StressState,
    velocity::VelocityState,
    prefactors::PseudoTransientSolver,
    ::MomentumBalance2D,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _membrane_stress_cached! => (
            stress.membrane_xx,
            stress.membrane_xy,
            stress.membrane_yy,
            prefactors.membrane_pre_aa,
            prefactors.membrane_pre_ab,
            velocity,
            mask,
        ),
    )
    return nothing
end

###############################################################
# Chmy-native, C-grid staggered Blatter-Pattyn membrane stress
###############################################################
#
# Per unit *volume*, not per unit area — no `H` anywhere, unlike the SSA/DIVA method above —
# with a genuine vertical-shear pair `σxz`/`σyz` the depth-integrated balance has none of:
#
#   σxx = 2µ(2u_x + v_y)   σxy = µ(u_y + v_x)   σxz = µ u_z
#
# and `σyy`/`σyz` the mirror image. Written into `stress.xx`/`xy`/`xz`/`yy`/`yz` directly:
# unlike SSA/DIVA's `membrane_xx`/`xy`/`yy`, no dedicated field is needed since BP's operand
# is already `AA3`/`AB3`/`ACXZ3`-shaped. `deviatoric_stress!` writes the true pointwise
# deviatoric stress into the same fields from `strainrate.xx`/`xy`/`xz` (a different
# combination, 2µε̇xx not 2µ(2ε̇xx+ε̇yy)) — the two are simply never live at once.
#
# `σxx`/`σyy` need no interpolation of `µ` (already at `aa`); `σxy` needs the ordinary
# one-way `hlerp` onto `ab`. `σxz`/`σyz` need `µ` staggered in *both* x and z onto
# `acx_ac`/`acy_ac`, which Chmy's `hlerp` handles with no special-casing — it interpolates
# every dimension on which `location(µ)` and the target differ.
#
# Masking is *not* uniform across the three sites, and the split is load-bearing. `σxy` at
# `ab` keeps the strict `node_fully_active` rule: `ab` is a genuine four-cell node, and a
# corner touching ice-free ground really does transmit no shear in this discretization.
# `σxz`/`σyz` do not, because `acx_ac`/`acy_ac` resolve to the *same* cell pair the unknown
# itself does — applying the strict rule there would delete BP's vertical operator along the
# whole margin ring. `_mu_acxz`/`_mu_acyz` carry that argument in full.

@kernel inbounds = true function _membrane_stress_staggered_bp!(
    sxx,
    sxy,
    sxz,
    syy,
    syz,
    μ,
    vel,
    mask,
    grid,
    O,
)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    Z = zero(eltype(sxx))

    if node_active(mask, NODE_AA, i, j)
        two_μ = 2 * μ[I...]
        dux = vel.x_dx[I...]
        dvy = vel.y_dy[I...]
        sxx[I...] = two_μ * (2 * dux + dvy)
        syy[I...] = two_μ * (dux + 2 * dvy)
    else
        sxx[I...] = Z
        syy[I...] = Z
    end

    sxy[I...] =
        node_fully_active(mask, NODE_AB, i, j) ?
        hlerp(μ, NODE_AB, grid, I...) * (vel.x_dy[I...] + vel.y_dx[I...]) : Z
    # `_mu_acxz`/`_mu_acyz` (`src/mechanics/pseudotransient.jl`), shared with the Gershgorin
    # bound rather than restated here: the two must apply the *same* `µ` to the same
    # interface or the bound stops describing the operator it bounds. They are one-sided at
    # the margin rather than strictly masked — see the note at their definition for why
    # `node_fully_active` is the wrong rule for this node class in particular.
    sxz[I...] = _mu_acxz(μ, mask, grid, I...) * vel.x_dz[I...]
    syz[I...] = _mu_acyz(μ, mask, grid, I...) * vel.y_dz[I...]
end

"""
$(TYPEDSIGNATURES)

Chmy-native, C-grid staggered Blatter-Pattyn membrane stress: writes `stress.xx`/`yy` (at
`aa`), `stress.xy` (at `ab`) and `stress.xz`/`yz` (at `acx`/`acy`, z-`Vertex`) from the 3D
viscosity `material.viscosity` and the column velocity gradients [`velocitygradients!`](@ref)
already wrote — `σxx = 2µ(2u_x+v_y)`, `σxy = µ(u_y+v_x)`, `σxz = µ u_z` (Robinson et al. 2022,
Eq. 1) and the `y` mirror image.

Per unit volume, unlike [`membranestress!(::MomentumBalance2D)`](@ref) — no thickness enters
anywhere, so no `topo` argument is taken. `µ` is interpolated onto `ab`/`acx_ac`/`acy_ac`
harmonically (`hlerp`), matching [`deviatoric_stress!`](@ref)'s stress-continuity argument;
`σxz`/`σyz` need `µ` staggered in *both* x and z, which `hlerp` handles with no extra code
(see the source note above).

Runs on `rt.grid`, the column grid — every operand is a genuine 3D field.

!!! note "`σxz`/`σyz` are one-sided at the margin, `σxy` is not"
    The `µ` behind `σxz`/`σyz` comes from `_mu_acxz`/`_mu_acyz`, shared with
    [`gershgorin_dt!`](@ref) so the bound and the operator cannot drift apart. Where only one
    of the two cells under an `acx`/`acy` face carries ice, `µ` is taken from that column
    alone rather than zeroed — see the note at their definition for why the strict rule that
    is right at `ab` is wrong here.

!!! warning "A zero viscosity gives `NaN`, not zero — same trap as `deviatoric_stress!`"
    `hlerp` averages reciprocals, so an unmasked ice-free neighbour produces `NaN` at `ab`
    rather than `0`; pass an [`IceMask`](@ref) once the material state has ice-free cells.
    The one-sided `acx_ac`/`acy_ac` path never inverts an ice-free cell's `µ` and so is safe
    either way.
"""
function membranestress!(
    stress::StressState,
    velocity::VelocityState,
    material::MechanicMaterialState,
    momentum::MomentumBalance3D,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    rt.launch(
        rt.arch,
        rt.grid,
        _membrane_stress_staggered_bp! => (
            stress.xx,
            stress.xy,
            stress.xz,
            stress.yy,
            stress.yz,
            material.viscosity,
            velocity,
            mask,
            rt.grid,
        ),
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level [`membranestress!`](@ref) for the Blatter-Pattyn momentum balance: writes
`mech.stress.xx`/`xy`/`xz`/`yy`/`yz` from `mech.material` and `mech.velocity` (which must
already carry the column velocity gradients, see [`velocitygradients!`](@ref)).
"""
membranestress!(
    mech::MechanicState,
    momentum::MomentumBalance3D,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
) = membranestress!(mech.stress, mech.velocity, mech.material, momentum, rt, mask)

###############################################################
# Chmy-native, C-grid staggered SSA/DIVA effective strain rate
###############################################################
#
# The true (not membrane-stress) second invariant ε̇_e of the *depth-averaged* velocity —
# Robinson et al. (2022), Eq. 12 — needed by `GlenViscosityContinuation`
# (`src/mechanics/solvers.jl`) to derive `viscosity_depthaveraged` from the current velocity
# iterate.
#
# Not reusable from `raw_strainrate_effective!`: that one reads `strainrate.xy` at `ab` as
# already written by `raw_strainrate!`, and works on column fields throughout.

@kernel inbounds = true function _effective_strainrate_ssa_staggered!(
    eff,
    velocity,
    mask,
    grid,
    O,
)
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
different quantity living in `strainrate.effective` (`AA3`), not the same field at two
resolutions.

Requires the depth-averaged velocity gradients to already be current — call
[`depthaverage_velocitygradients!`](@ref) first, exactly as [`membranestress!`](@ref) does.
"""
function effective_strainrate_ssa!(
    strainrate::StrainRateState,
    velocity::VelocityState,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _effective_strainrate_ssa_staggered! =>
            (strainrate.effective_depthaveraged, velocity, mask, rt.grid2d),
    )
    return nothing
end

###############################################################
# Chmy-native DIVA effective strain rate (per layer)
###############################################################
#
# See the docstring below for the equations; implementation notes not there:
#
# `τ_b` lives on the velocity faces (`acx`/`acy`) and `ε̇_e` at `aa`, so each component is
# `lerp`ed onto `aa` — a 2D interpolation reused down the whole column rather than a 3D one
# per layer. The horizontal invariant is likewise recomputed per layer instead of cached in
# a 2D scratch field: a handful of flops against a field allocation, and it keeps the kernel
# a pure function of state. Consequently DIVA never reads `velocity.x_dz`.

@kernel inbounds = true function _effective_strainrate_diva!(
    eff,
    velocity,
    μ,
    τbx,
    τby,
    mask,
    grid,
    grid2d,
    O,
)
    I = @index(Global, NTuple)
    I = I + O
    i, j, k = I
    T = eltype(eff)
    Z = zero(T)

    if node_active(mask, NODE_AA, i, j)
        # Horizontal part (Eq. 12): depth-independent, from the depth-averaged gradients.
        dxx = velocity.depthaverage_x_dx[i, j, 1]
        dyy = velocity.depthaverage_y_dy[i, j, 1]
        dxy = lerp(velocity.depthaverage_x_dy, NODE_AA, grid2d, i, j, 1)
        dyx = lerp(velocity.depthaverage_y_dx, NODE_AA, grid2d, i, j, 1)
        horizontal = dxx^2 + dyy^2 + dxx * dyy + ((dxy + dyx) / 2)^2

        # Vertical shear (Eq. 21), with H already cancelled — see the note above.
        w = one(T) - zcenter(grid, k)               # (s - z)/H at the layer midpoint
        μk = μ[i, j, k]
        uz = μk > Z ? lerp(τbx, NODE_AA, grid2d, i, j, 1) * w / μk : Z
        vz = μk > Z ? lerp(τby, NODE_AA, grid2d, i, j, 1) * w / μk : Z

        eff[I...] = sqrt(horizontal + (uz^2 + vz^2) / 4)
    else
        eff[I...] = Z
    end
end

"""
$(TYPEDSIGNATURES)

DIVA's per-layer effective strain rate (Robinson et al. 2022, Eq. 13), written to
`strainrate.effective` (`AA3`) on the column grid:

```math
\\dot\\varepsilon_e^2 = \\bar u_x^2 + \\bar v_y^2 + \\bar u_x \\bar v_y
    + \\tfrac14(\\bar u_y + \\bar v_x)^2 + \\tfrac14 u_z^2 + \\tfrac14 v_z^2
```

The horizontal terms are the SSA invariant of the *depth-averaged* velocity, identical in
every layer; the depth dependence enters only through the vertical shear, which is
**diagnosed from the basal stress** via Eq. (21) rather than differentiated from a velocity
profile:

```math
u_z(z) = \\frac{\\tau_{b,x}\\,(s-z)}{\\eta(z)\\,H} = \\frac{\\tau_{b,x}\\,(1-\\zeta)}{\\mu(z)}
```

(the thickness cancels on the sigma axis). This is the only form available inside the
pseudo-transient loop, where no 3D velocity exists yet — `τ_b` and `µ` are the previous
iterate's, as the paper prescribes.

Contrast [`effective_strainrate_ssa!`](@ref), which computes Eq. (12) — the same expression
with the shear terms dropped — into `strainrate.effective_depthaveraged` (`AA2`).

Requires the depth-averaged velocity gradients ([`depthaverage_velocitygradients!`](@ref))
and the basal stress ([`basalstress!`](@ref)) to be current, and `material.viscosity` to
hold the previous iterate's `µ(z)`.
"""
function effective_strainrate_diva!(
    strainrate::StrainRateState,
    velocity::VelocityState,
    material::MechanicMaterialState,
    stress::StressState,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    rt.launch(
        rt.arch,
        rt.grid,
        _effective_strainrate_diva! => (
            strainrate.effective,
            velocity,
            material.viscosity,
            stress.base_x,
            stress.base_y,
            mask,
            rt.grid,
            rt.grid2d,
        ),
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level [`effective_strainrate_diva!`](@ref).
"""
effective_strainrate_diva!(
    mech::MechanicState,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
) = effective_strainrate_diva!(
    mech.strainrate,
    mech.velocity,
    mech.material,
    mech.stress,
    rt,
    mask,
)

###############################################################
# Chmy-native Blatter-Pattyn effective strain rate (per layer)
###############################################################
#
# Same invariant as DIVA (Eq. 13), but every term is genuinely 3D: BP already carries a
# resolved column velocity, so `u_z`/`v_z` are read straight from `velocity.x_dz`/`y_dz`
# rather than diagnosed from `τ_b`. No `µ`, no basal stress and no `grid2d` argument are
# needed here — the sole difference from `_effective_strainrate_diva!`.
#
# `x_dx`/`y_dy` are already at `aa`; `x_dy`/`y_dx` at `ab` and `x_dz`/`y_dz` at
# `acx_ac`/`acy_ac` each need one `lerp` onto `aa` (the latter a two-way x-and-z stagger,
# which Chmy's `lerp` handles with no special-casing). Arithmetic, not harmonic: this
# interpolates a strain rate, not a viscosity, so there's no NaN trap either —
# `_velocity_gradients!` already zeroes an inactive neighbour rather than leaving a
# placeholder to invert.

@kernel inbounds = true function _effective_strainrate_bp!(eff, velocity, mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    if node_active(mask, NODE_AA, i, j)
        dxx = velocity.x_dx[I...]
        dyy = velocity.y_dy[I...]
        dxy = lerp(velocity.x_dy, NODE_AA, grid, I...)
        dyx = lerp(velocity.y_dx, NODE_AA, grid, I...)
        uz = lerp(velocity.x_dz, NODE_AA, grid, I...)
        vz = lerp(velocity.y_dz, NODE_AA, grid, I...)
        eff[I...] =
            sqrt(dxx^2 + dyy^2 + dxx * dyy + ((dxy + dyx) / 2)^2 + (uz^2 + vz^2) / 4)
    else
        eff[I...] = zero(eltype(eff))
    end
end

"""
$(TYPEDSIGNATURES)

Chmy-native Blatter-Pattyn effective strain rate at `aa` (Robinson et al. 2022, Eq. 3):

```math
\\dot\\varepsilon_e^2 = u_x^2 + v_y^2 + u_x v_y + \\tfrac14(u_y + v_x)^2
    + \\tfrac14 u_z^2 + \\tfrac14 v_z^2
```

Writes `strainrate.effective` (`AA3`), the same field [`effective_strainrate_diva!`](@ref)
writes — the two momentum balances never run in the same solve, so sharing the field costs
nothing. Unlike DIVA's Eq. (13), every term (including `u_z`/`v_z`) comes from the real
velocity gradients [`velocitygradients!`](@ref) already wrote, since BP resolves the column
velocity directly rather than reconstructing it afterwards.

Requires `velocitygradients!` to have already run (all nine gradient fields current).
"""
function effective_strainrate_bp!(
    strainrate::StrainRateState,
    velocity::VelocityState,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    rt.launch(
        rt.arch,
        rt.grid,
        _effective_strainrate_bp! => (strainrate.effective, velocity, mask, rt.grid),
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level [`effective_strainrate_bp!`](@ref).
"""
effective_strainrate_bp!(
    mech::MechanicState,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
) = effective_strainrate_bp!(mech.strainrate, mech.velocity, rt, mask)

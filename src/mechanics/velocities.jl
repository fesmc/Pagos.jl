"""
$(TYPEDSIGNATURES)

Compute the basal velocity `vb` from the depth-averaged velocity `v`.
"""
function basal_velocity_from_depthavg_velocity!(vb, v, beta, F2)
    @. vb = v / (1 + beta * F2)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Compute the basal velocity `vb` from the surface velocity `vs`.
"""
function basal_velocity_from_surface_velocity!(vb, vs, beta, F1)
    @. vb = vs / (1 + beta * F1)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Compute the surface velocity `vs` from the basal velocity `vb`.
"""
function surface_velocity!(vs, vb, beta, F1)
    @. vs = vb * (1 + beta * F1)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Compute the depth-averaged velocity `v` from the basal velocity `vb`.
"""
function depthavg_velocity!(v, vb, beta, F2)
    @. v = vb * (1 + beta * F2)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Compute the 3D velocity field. `F1` is used as a work array for the running viscosity
integral; it is zeroed on entry. `sigma` must be a host-side vector of the sigma-coordinate
**upper interface** of each layer — `ζ_ac[2:end]`, i.e. `nz` ascending values ending at
`1`, as produced by `exponential_vertical_layers` — accessible by index on the CPU
during the layer loop. *Not* the layer midpoints `ζ_aa`: the layer thickness is recovered
here as `sigma[l] - sigma[l-1]` (with `sigma[0] ≡ 0`), which is a half-layer short for the
first layer if midpoints are passed, and the midpoint is then reconstructed as
`sigma[l] - dsigma/2`.

Because the running integral is advanced across the *whole* of layer `l` before the layer
velocity is written, `v[:, :, l]` is the velocity at that layer's **upper interface**, not
at its midpoint — which is why `v[:, :, end]` is the surface velocity.
"""
function velocities3D!(F1, vx3D, vy3D, vb_x, vb_y, mu, beta, H, sigma)
    F1 .= zero(eltype(F1))
    for l in eachindex(sigma)
        aggregate_viscosity_integral!(F1, mu, H, 1, sigma, l)
        layer_velocity!(vx3D, l, vb_x, beta, F1)
        layer_velocity!(vy3D, l, vb_y, beta, F1)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Aggregate the viscosity integral over all sigma layers into `Fm`. `sigma` must be a
host-side vector of layer **upper interfaces** (see [`velocities3D!`](@ref)).
"""
function aggregated_viscosity_integral!(Fm, mu, H, m, sigma)
    Fm .= zero(eltype(Fm))
    for l in eachindex(sigma)
        aggregate_viscosity_integral!(Fm, mu, H, m, sigma, l)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Accumulate the `l`-th layer's contribution to the viscosity integral `Fm` using a
midpoint Riemann sum in sigma coordinates — Robinson et al. (2022), Eq. 15. `sigma` must be
a host-side vector of layer **upper interfaces** (see [`velocities3D!`](@ref)); the layer
thickness `dsigma` and the midpoint value of `(s - z)/H` are derived from it here.
GPU-compatible.
"""
function aggregate_viscosity_integral!(Fm, mu, H, m, sigma, l)
    dsigma = l == 1 ? sigma[l] : sigma[l] - sigma[l-1]
    s_minus_z_over_H = 1 - sigma[l] + dsigma / 2
    backend = get_backend(Fm)
    kernel! = _aggregate_viscosity_integral!(backend)
    kernel!(Fm, mu, H, m, s_minus_z_over_H, dsigma, l; ndrange = size(Fm))
    return nothing
end

"""
$(TYPEDSIGNATURES)

Write the 3D velocity at layer `l` into `v` given the basal velocity `vb` and the
viscosity integral `F1` accumulated up to that layer. GPU-compatible.
"""
function layer_velocity!(v, l, vb, beta, F1)
    backend = get_backend(v)
    kernel! = _layer_velocity!(backend)
    kernel!(v, vb, beta, F1, l; ndrange = size(vb))
    return nothing
end

@kernel function _aggregate_viscosity_integral!(
    Fm,
    mu,
    H,
    m,
    s_minus_z_over_H,
    dsigma,
    l,
)
    i, j = @index(Global, NTuple)
    @inbounds Fm[i, j] += (s_minus_z_over_H ^ m * dsigma * H[i, j]) / mu[i, j, l]
end

@kernel function _layer_velocity!(v, vb, beta, F1, l)
    i, j = @index(Global, NTuple)
    @inbounds v[i, j, l] = vb[i, j] * (1 + beta[i, j] * F1[i, j])
end

###############################################################
# Chmy-native, C-grid staggered DIVA viscosity integrals
###############################################################
#
# The σ values come from the grid, not the caller: `zcenter`/`Δz` are Chmy's own cell
# midpoint and interface-to-interface thickness, so the midpoint rule is exact-by-
# construction on any layering — unlike the collocated path below, whose `sigma` argument
# must be layer upper interfaces and silently mis-weights the basal layer if midpoints are
# passed instead (a trap that can't be reached from here). See the docstring for the maths.
#
# `F_1` and `F_2` are fused into one pass since they share `Δζ_k`, the thickness scaling and
# the expensive part, the reciprocal viscosity — the only two orders DIVA needs (`F_1` for
# the surface velocity Eq. 17, `F_2` for the depth-averaged velocity/`β_eff`, Eqs. 18-19).
#
# One thread per column: the launch is 2D, with the vertical sum serial inside the kernel.

@kernel inbounds = true function _viscosity_integrals!(F1, F2, μ, H, nz, mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    T = eltype(F1)
    Z = zero(T)
    # `H`, `F1`, `F2` are `grid2d` fields read at the sweep's own index (no staggering);
    # `μ` is a column field, so its `k` is this kernel's own loop variable, not the sweep's.
    Hij = H[I...]

    if node_active(mask, NODE_AA, i, j) && Hij > Z
        f1 = Z
        f2 = Z
        for k = 1:nz
            w = one(T) - zcenter(grid, k)              # (s - z)/H at the layer midpoint
            r = Δz(grid, Center(), i, j, k) * Hij / μ[i, j, k]   # dz / µ
            f1 += w * r
            f2 += w * w * r
        end
        F1[I...] = f1
        F2[I...] = f2
    else
        F1[I...] = Z
        F2[I...] = Z
    end
end

"""
$(TYPEDSIGNATURES)

Chmy-native, C-grid staggered DIVA viscosity integrals (Robinson et al. 2022, Eq. 15):
fill `F1` and `F2` at `aa` on `rt.grid2d` with

```math
F_m = \\int_b^s \\frac{1}{\\mu} \\left(\\frac{s - z}{H}\\right)^m \\mathrm{d}z
```

for `m = 1, 2`, from the 3D viscosity `viscosity` (at `aa` on the column grid `rt.grid`)
and the thickness `H` (at `aa` on `rt.grid2d`). `F1` sets the surface velocity
`u_s = u_b (1 + βF₁)` (Eq. 17) and `F2` the depth-averaged velocity `ū = u_b (1 + βF₂)`
(Eq. 18) and hence the DIVA effective friction `β_eff = β/(1 + βF₂)` (Eq. 19).

Distinguished from the collocated [`aggregated_viscosity_integral!`](@ref) by taking a
[`Runtime`](@ref) — and by taking **no `sigma` argument at all**: the layer midpoints and
thicknesses are read from the grid's own sigma axis, so there is no interface-vs-midpoint
convention for a caller to get wrong (see the source note above). Both integrals are
computed in a single pass over the column.

Ice-free columns (`H == 0`, or masked out) get `0`, not a division by zero.

# Accuracy

The midpoint rule is exact for integrands linear in `σ`, so under a vertically uniform
viscosity `F₁ = H/(2µ)` comes out exact to roundoff **on any layering**, however stretched.
`F₂ = H/(3µ)` is second-order in the layer thickness, as is `F₁` once `µ` varies with depth.

!!! warning "Needs a real column: `nz == 1` cannot resolve `F₂`"
    On a depth-integrated grid (`nz == 1`, `rt.grid === rt.grid2d`) the single layer puts
    the quadrature point at `σ = 1/2`, giving `F₂ = H/(4µ)` against the true `H/(3µ)` — 25%
    low, and no amount of masking or halo care changes that. `F₁` stays exact. DIVA is a
    column approximation; this function is only meaningful on a grid that has a column.

!!! warning "A zero viscosity gives `Inf`, not zero"
    The integrand is `1/µ`, so an unmasked ice-free layer produces `Inf` rather than `0` —
    mechanically the same trap as `hlerp`'s `NaN` in [`strainrate!`](@ref). The `H > 0`
    guard catches the ordinary ice-free column; pass an [`IceMask`](@ref) once the material
    state carries ice-free cells with a placeholder viscosity.

!!! warning "`viscosity` and `H` need their halo filled by the caller"
    Both are read at the sweep's own index, which the `Launcher` extends one ring beyond
    the interior — the standard Chmy-native contract (`setdata!` fills the interior only;
    see [`pseudo_transient!`](@ref)'s note). Unlike the stencil kernels there is no
    interpolation here, so an unfilled halo costs only the halo ring of `F1`/`F2`, not the
    interior.
"""
function viscosity_integrals!(
    F1,
    F2,
    viscosity,
    H,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    nz = size(rt.grid, Center())[3]
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _viscosity_integrals! => (F1, F2, viscosity, H, nz, mask, rt.grid),
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level [`viscosity_integrals!`](@ref): fill `F1`/`F2` from `mech.material.viscosity`
(the 3D viscosity, *not* the depth-averaged one) and `mech.topography.thickness`.

`F1`/`F2` take no home in [`MechanicState`](@ref) yet, so they stay explicit arguments.
"""
viscosity_integrals!(
    F1,
    F2,
    mech::MechanicState,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
) = viscosity_integrals!(
    F1,
    F2,
    mech.material.viscosity,
    mech.topography.thickness,
    rt,
    mask,
)

###############################################################
# Column → depth-average reduction
###############################################################
#
# `f̄ = (1/H) ∫_b^s f dz`. On the sigma axis `dz = H dζ`, so the thickness cancels and the
# depth average is a pure `ζ` quadrature: `f̄ = ∫_0^1 f dζ ≈ Σ_k f[i,j,k] Δζ_k`, with the
# weights summing to 1 by construction. That is why this takes no `H` — an average, unlike
# the `F_m` *integrals* above, is thickness-independent.

@kernel inbounds = true function _depthaverage!(out, f, nz, mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    T = eltype(out)
    acc = zero(T)
    if node_active(mask, NODE_AA, i, j)
        for k = 1:nz
            acc += f[i, j, k] * Δz(grid, Center(), i, j, k)
        end
    end
    out[I...] = acc
end

"""
$(TYPEDSIGNATURES)

Depth-average the column field `f` (at `aa` on `rt.grid`) into `out` (at `aa` on
`rt.grid2d`): `f̄ = (1/H)∫_b^s f dz = ∫_0^1 f dζ`, the layer-thickness-weighted mean over
the sigma axis.

Takes no thickness argument: `dz = H dζ` makes the `H` cancel in an *average* (contrast
[`viscosity_integrals!`](@ref), where the integral keeps a factor of `H`). The weights
`Δζ_k` come from the grid's own sigma axis and sum to exactly 1, so a vertically uniform
field reproduces its own value to roundoff on any layering.

DIVA's use is `µ̄`, the depth-averaged viscosity the membrane stress needs (Robinson et al.
2022, Eq. 14) — but the reduction is generic.
"""
function depthaverage!(out, f, rt::Runtime, mask::AbstractIceMask = NoMask())
    nz = size(rt.grid, Center())[3]
    rt.launch2d(rt.arch, rt.grid2d, _depthaverage! => (out, f, nz, mask, rt.grid))
    return nothing
end

###############################################################
# DIVA effective basal friction
###############################################################
#
# Implemented as the algebraically identical reciprocal form `β_eff = 1/(1/β + F₂)` rather
# than the paper's `β/(1+βF₂)` (Eq. 19) — see the docstring below for why that covers both
# the frozen-bed and free-slip limits with no branch.

@kernel inbounds = true function _beta_eff_diva!(β_eff, β, F2, mask, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    T = eltype(β_eff)
    if node_active(mask, NODE_AA, i, j)
        # `inv(0) == Inf` and `inv(Inf) == 0` are exactly the two limits wanted here.
        β_eff[I...] = inv(inv(β[I...]) + F2[I...])
    else
        β_eff[I...] = zero(T)
    end
end

"""
$(TYPEDSIGNATURES)

DIVA's effective basal friction coefficient (Robinson et al. 2022, Eq. 19):

```math
\\beta_{\\mathrm{eff}} = \\frac{\\beta}{1 + \\beta F_2}
    \\;=\\; \\frac{1}{1/\\beta + F_2}
```

written into `beta_eff` from the friction coefficient `beta` and the viscosity integral
`F2` (see [`viscosity_integrals!`](@ref)). It is what makes the basal stress expressible
against the depth-averaged velocity, `τ_b = β_eff·ū`, which is the whole reason DIVA can be
solved as a depth-integrated balance.

The reciprocal form is used deliberately: it reproduces the paper's frozen-bed limit
Eq. (20), `β_eff = 1/F₂` as `β → ∞`, and the free-slip case `β_eff = 0` at `β = 0`, both
through ordinary IEEE arithmetic rather than a branch on a tunable "large β" threshold.

!!! note "SSA is the `F₂ = 0` limit"
    With `F₂ = 0` this returns `β` unchanged, i.e. `τ_b = βū` — the SSA friction law. So a
    zeroed `F₂` degrades to SSA rather than to nonsense, which is the sane failure mode if
    the integrals were never computed.
"""
function beta_eff_diva!(
    beta_eff,
    beta,
    F2,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    rt.launch2d(rt.arch, rt.grid2d, _beta_eff_diva! => (beta_eff, beta, F2, mask))
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level [`beta_eff_diva!`](@ref): writes `mech.friction.beta_eff` from
`mech.friction.beta` and `mech.material.viscosity_integral_2`.
"""
beta_eff_diva!(mech::MechanicState, rt::Runtime, mask::AbstractIceMask = NoMask()) =
    beta_eff_diva!(
        mech.friction.beta_eff,
        mech.friction.beta,
        mech.material.viscosity_integral_2,
        rt,
        mask,
    )

###############################################################
# Chmy-native, C-grid staggered 3D velocity reconstruction
###############################################################
#
# Eq. 16 (with `H` included, unlike `effective_strainrate_diva!`'s Eq. 21 where it cancels):
#
#   u(z) = u_b + (βu_b/H) ∫_b^z (s-z')/µ(z') dz'
#
# On the sigma axis, writing F₁(ζ) for the *partial* (running) form of the Eq. 15 integral:
#
#   u(ζ) = u_b·(1 + β·F₁(ζ)),   F₁(ζ) = (H/µ)∫_0^ζ (1-σ) dσ = (H/µ)·(ζ - ζ²/2)
#
# the same relation as Eq. 17 (`u_s = u_b(1+βF₁)`), stopped at an interior ζ instead of the
# surface. `u_b` is Eq. 18 inverted, `u_b = ū/(1 + βF₂) = ū·β_eff/β`. `β` (bare), not
# `β_eff`, is what Eq. 16/18 are written in terms of.
#
# Evaluated at the layer's own z-`Center`, not its upper interface (contrast the collocated
# `velocities3D!` below, whose running sum crosses a whole layer before writing and lands on
# z-`Vertex`). `g(σ) = σ - σ²/2`, the exact antiderivative of `(1-σ)`, makes the per-layer
# partial and full contributions exact rather than a quadrature approximation, under the
# same piecewise-constant-per-layer `µ` every `F_m` integral in this file assumes — a naive
# "half of the full layer's contribution" would *not* be exact here since `g` is quadratic.
#
# One thread per column, matching `viscosity_integrals!`: the running sum is inherently
# serial, so there's nothing to parallelize within a column.

@kernel inbounds = true function _velocities3D_diva!(
    vx,
    vy,
    ubar_x,
    ubar_y,
    β,
    F2,
    μ,
    H,
    nz,
    mask,
    grid,
    O,
)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    T = eltype(vx)
    Z = zero(T)
    Hij = H[i, j, 1]

    if node_active(mask, NODE_AA, i, j) && Hij > Z
        βij = β[i, j, 1]
        ubx = ubar_x[i, j, 1] / (one(T) + βij * F2[i, j, 1])   # Eq. 18 inverted
        uby = ubar_y[i, j, 1] / (one(T) + βij * F2[i, j, 1])
        g(σ) = σ - σ^2 / 2
        running = Z                       # F₁ accumulated up to the *previous* interface
        prev = zvertex(grid, 1)
        for k = 1:nz
            mid = zcenter(grid, k)
            next = zvertex(grid, k + 1)
            invμ = inv(μ[i, j, k])
            partial = (g(mid) - g(prev)) * Hij * invμ
            full = (g(next) - g(prev)) * Hij * invμ
            vx[i, j, k] = ubx * (one(T) + βij * (running + partial))
            vy[i, j, k] = uby * (one(T) + βij * (running + partial))
            running += full
            prev = next
        end
    else
        for k = 1:nz
            vx[i, j, k] = Z
            vy[i, j, k] = Z
        end
    end
end

# Surface velocity, Eq. 17, `u_s = u_b(1+βF₁)`: a plain 2D kernel reusing the already
# computed full-column `F₁`/`F₂` (`material.viscosity_integral_1`/`_2`) rather than reading
# off the top of `_velocities3D_diva!`'s profile, which sits at the last layer's midpoint,
# not at ζ = 1. No column loop needed.
@kernel inbounds = true function _surfacevelocity_diva!(
    vsx,
    vsy,
    ubar_x,
    ubar_y,
    β,
    F1,
    F2,
    mask,
    O,
)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    T = eltype(vsx)
    if node_active(mask, NODE_AA, i, j)
        βij = β[I...]
        ubx = ubar_x[I...] / (one(T) + βij * F2[I...])
        uby = ubar_y[I...] / (one(T) + βij * F2[I...])
        vsx[I...] = ubx * (one(T) + βij * F1[I...])
        vsy[I...] = uby * (one(T) + βij * F1[I...])
    else
        vsx[I...] = zero(T)
        vsy[I...] = zero(T)
    end
end

"""
$(TYPEDSIGNATURES)

Chmy-native, C-grid staggered 3D velocity reconstruction under [`DIVAMomentumBalance`](@ref)
(Robinson et al. 2022, Eq. 16): writes `velocity.x`/`y` (at z-`Center`, one value per layer)
from the depth-averaged solve `velocity.depthaverage_x`/`y`, the friction coefficient
`friction.beta` and the viscosity integrals `material.viscosity_integral_1`/`_2` — all of
which [`diva_update!`](@ref) is responsible for having populated. Also writes
`velocity.surface_x`/`y` (Eq. 17).

A **diagnostic**, not part of the momentum solve: call this *after* [`pseudo_transient!`](@ref)
converges, on whatever state it left. `pseudo_transient!` never calls it itself.

Distinguished from the collocated [`velocities3D!`](@ref) by taking a [`Runtime`](@ref) and
a [`DIVAMomentumBalance`](@ref) — and, more substantively, by evaluating each layer at its
own z-`Center` (matching `velocity.x`'s own node class) rather than at the layer's upper
z-`Vertex` interface, which is what the collocated method's running-sum-then-write order
produces. The two are not comparable point for point on a shared grid; see the source note
above.

!!! note "Requires `nz > 1`"
    Trivially true whenever [`DIVAMomentumBalance`](@ref) itself is usable at all — see
    [`pseudo_transient!`](@ref)'s `nz == 1` guard.
"""
function velocities3D!(
    mech::MechanicState,
    rt::Runtime,
    ::DIVAMomentumBalance,
    mask::AbstractIceMask = NoMask(),
)
    (; velocity, friction, material, topography) = mech
    nz = size(rt.grid, Center())[3]
    # `rt.launch2d`, not `rt.launch`: the vertical dependency is a serial running sum, so
    # this wants one thread per *column* (as `viscosity_integrals!` already established),
    # not one thread per `(i,j,k)` triple — the latter would have every thread in a column
    # redundantly re-run the same full-column loop.
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _velocities3D_diva! => (
            velocity.x,
            velocity.y,
            velocity.depthaverage_x,
            velocity.depthaverage_y,
            friction.beta,
            material.viscosity_integral_2,
            material.viscosity,
            topography.thickness,
            nz,
            mask,
            rt.grid,
        ),
    )
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _surfacevelocity_diva! => (
            velocity.surface_x,
            velocity.surface_y,
            velocity.depthaverage_x,
            velocity.depthaverage_y,
            friction.beta,
            material.viscosity_integral_1,
            material.viscosity_integral_2,
            mask,
        ),
    )
    return nothing
end

# SSA is plug flow: u(z) = ū at every layer, no shear at all, so the "reconstruction" is a
# broadcast, not a quadrature. Included so that a caller reading `velocity.x`/`y`/
# `surface_x`/`y` after a solve never has to know which momentum balance produced them.
@kernel inbounds = true function _velocities3D_ssa!(vx, vy, ubar_x, ubar_y, nz, mask, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    T = eltype(vx)
    ux = node_active(mask, NODE_AA, i, j) ? ubar_x[i, j, 1] : zero(T)
    uy = node_active(mask, NODE_AA, i, j) ? ubar_y[i, j, 1] : zero(T)
    for k = 1:nz
        vx[i, j, k] = ux
        vy[i, j, k] = uy
    end
end

"""
$(TYPEDSIGNATURES)

Chmy-native 3D velocity "reconstruction" under [`SSAMomentumBalance`](@ref): plug flow, so
`velocity.x`/`y` and `velocity.surface_x`/`y` are simply set to the depth-averaged solve at
every layer/at the surface — no shear, no viscosity integrals read. Exists so that
`velocity.x`/`y`/`surface_x`/`y` are always populated after a solve, regardless of which
[`MomentumBalance2D`](@ref) balance produced it.
"""
function velocities3D!(
    mech::MechanicState,
    rt::Runtime,
    ::SSAMomentumBalance,
    mask::AbstractIceMask = NoMask(),
)
    (; velocity) = mech
    nz = size(rt.grid, Center())[3]
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _velocities3D_ssa! => (
            velocity.x,
            velocity.y,
            velocity.depthaverage_x,
            velocity.depthaverage_y,
            nz,
            mask,
        ),
    )
    copyto!(asarray(velocity.surface_x), asarray(velocity.depthaverage_x))
    copyto!(asarray(velocity.surface_y), asarray(velocity.depthaverage_y))
    return nothing
end

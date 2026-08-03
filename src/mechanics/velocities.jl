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
    dsigma = l == 1 ? sigma[l] : sigma[l] - sigma[l - 1]
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

@kernel function _aggregate_viscosity_integral!(Fm, mu, H, m, s_minus_z_over_H, dsigma, l)
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
# Robinson et al. (2022), Eq. 15 — the generalized integrals DIVA writes its basal stress
# and its vertical velocity profile in terms of:
#
#   F_m ≡ ∫_b^s (1/µ) ((s - z)/H)^m dz
#
# On the terrain-following sigma axis σ = (z - b)/H ∈ [0, 1] (bed at 0, surface at 1, see
# `AbstractSigmaTransform`), `(s - z)/H = 1 - σ` and `dz = H dσ`, so the whole integral is
# a pure σ quadrature scaled by the local thickness:
#
#   F_m = H ∫_0^1 (1 - σ)^m / µ(σ) dσ  ≈  H Σ_k (1 - ζ_aa[k])^m Δζ_k / µ[i, j, k]
#
# the midpoint rule over the grid's own layers, exactly as the collocated
# `aggregate_viscosity_integral!` above.
#
# **The σ values come from the grid, not from the caller.** `zcenter(grid, k)` *is* Chmy's
# cell midpoint and `Δz(grid, Center(), ...)` *is* the interface-to-interface thickness
# (`_sigma_axis`'s ghost interfaces make both total at the bed and surface), so the
# midpoint rule is exact-by-construction on any layering rather than by agreement with a
# convention. This is deliberately unlike the collocated path, whose `sigma` argument must
# be layer *upper interfaces* and silently mis-weights the basal layer if midpoints are
# passed instead — a trap that cannot be reached from here.
#
# `F_1` and `F_2` are fused into one pass: they share `Δζ_k`, the thickness scaling and —
# the expensive part — the reciprocal viscosity, so the second integral is nearly free.
# Those are also the only two orders DIVA ever needs (`F_1` for the surface velocity,
# Eq. 17; `F_2` for the depth-averaged velocity and hence `β_eff`, Eqs. 18–19).
#
# One thread per column: the launch is 2D (the output is), with the vertical sum serial
# inside the kernel. No atomics, no `nz` intermediate storage.

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
        for k in 1:nz
            w  = one(T) - zcenter(grid, k)              # (s - z)/H at the layer midpoint
            r  = Δz(grid, Center(), i, j, k) * Hij / μ[i, j, k]   # dz / µ
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
function viscosity_integrals!(F1, F2, viscosity, H, rt::Runtime,
                              mask::AbstractIceMask = NoMask())
    nz = size(rt.grid, Center())[3]
    rt.launch2d(rt.arch, rt.grid2d,
                _viscosity_integrals! => (F1, F2, viscosity, H, nz, mask, rt.grid))
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level [`viscosity_integrals!`](@ref): fill `F1`/`F2` from `mech.material.viscosity`
(the 3D viscosity, *not* the depth-averaged one) and `mech.topography.thickness`.

`F1`/`F2` take no home in [`MechanicState`](@ref) yet — the DIVA path that consumes them is
still being ported (`roadmaps/chmy.md`, Phase 3) — so they stay explicit arguments.
"""
viscosity_integrals!(F1, F2, mech::MechanicState, rt::Runtime,
                     mask::AbstractIceMask = NoMask()) =
    viscosity_integrals!(F1, F2, mech.material.viscosity, mech.topography.thickness,
                         rt, mask)

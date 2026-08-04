###############################################################
# Pseudo-transient momentum solver (Sandip et al., 2024)
###############################################################

"""
$(TYPEDSIGNATURES)

Relax the velocity field: `v = v_old + theta_v * dotvel * dtau`.
"""
function pseudo_vel!(v, v_old, dotvel, dtau, theta_v)
    @. v = v_old + theta_v * dotvel * dtau
    return nothing
end

###############################################################
# C-grid staggered pseudo-transient solver
###############################################################
#
# Two facts govern everything below and neither lives inside a kernel. `Chmy.Field` is an
# `AbstractArray` with no `BroadcastStyle`, so `copyto!`, the `v = v_old + θ·dv·dτ`
# relaxation and the max-norm reductions would fall back to scalar `getindex` — right on
# CPU, catastrophic on GPU — and are therefore all routed through `asarray`. And the
# membrane stress reads velocity gradients one ring beyond the interior at `ab`, so the
# velocity halo must be refreshed every iteration, not once at the end; `Neumann(0)` is the
# placeholder until per-equation boundary conditions land (`roadmaps/chmy.md` Phase 4).

"""
$(TYPEDSIGNATURES)

Local (per-grid-point) pseudo-transient time step (Sandip et al. 2024, Eq. 7): writes
`dtau_x`, `dtau_y` from the depth-averaged viscosity `mu` (`aa`), staggered onto the
velocity's own node class (`acx`/`acy`) by `lerp` — the same arithmetic-mean choice
[`drivingstress!`](@ref) and [`dotvel!`](@ref) make for the ice thickness `H`, and folding
the safety factor `dtau_scaling` in up front.

Per grid point rather than a `maximum(mu)` reduction, so a single stiff cell no longer
throttles `Δτ` everywhere else and no global reduction is needed. Only `mu`'s interior is
read (`lerp` needs its halo only at the two domain-boundary `ab`-style corners the membrane
stress itself needs — see [`pseudo_transient!`](@ref)'s halo warning), so this call has the
same halo requirement as the rest of the solve.

!!! note "`ρ` here is a numerical pseudo-density, not ice"
    Callers pass `c.density_ice`, but this `Δτ` governs an artificial dynamic relaxation,
    not physical inertia — `ρ` only sets a scale that `dtau_scaling` is free to absorb. So
    it is *not* something the `(m, yr, Pa)` convention (see [`Constants`](@ref)) obliges
    you to rescale: the viscosities this has always been tuned against were already
    `Pa yr`, and the value works as-is. Retune via `dtau_scaling` if a problem needs it.
"""
function pseudo_dt!(dtau_x, dtau_y, ρ, dx, dy, mu, muB, ndim, dtau_scaling, rt::Runtime)
    scaling = convert(eltype(dtau_x), dtau_scaling * ρ * dx * dy / (4 * (1 + muB) * ndim))
    rt.launch2d(rt.arch, rt.grid2d,
              _pseudo_dt_staggered! => (dtau_x, dtau_y, mu, scaling, rt.grid2d))
    return nothing
end

@kernel inbounds = true function _pseudo_dt_staggered!(dtau_x, dtau_y, mu, scaling, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    dtau_x[I...] = scaling / lerp(mu, NODE_ACX, grid, I...)
    dtau_y[I...] = scaling / lerp(mu, NODE_ACY, grid, I...)
end

###############################################################
# Gershgorin pseudo-time step
###############################################################
#
# `ViscosityPseudoTimeStep` (Sandip Eq. 7) bounds the *membrane* part of the operator and
# nothing else. That is fine on the uniform-slab tests, where β is small and there is no ice
# margin, and wrong in both directions on real geometry:
#
#  - it omits the basal-drag term β/(ρH) from the bound, which under grounded Antarctic ice
#    (β_eff up to ~5e13 Pa s m⁻¹) is the *dominant* eigenvalue — a solve that recomputes
#    friction from the velocity (`ActiveFrictionUpdate`) diverges within ~20 iterations;
#  - it `lerp`s the viscosity straight across the ice margin, so an ice-free cell's
#    placeholder viscosity throttles Δτ on exactly the faces that carry the calving-front
#    forcing.
#
# The bound below is Gershgorin's: λ_max ≤ max row sum of |coefficients| of the mass-scaled
# residual operator, and explicit stability is Δτ ≤ 2/λ_max. The coefficients are read off
# the discretisation the solver actually runs — `_membrane_stress_staggered!`'s `ηH` at `aa`
# and `hlerp(η)·lerp(H)` at `ab`, `_basalstress_staggered!`'s `lerp(β)`, and
# `_dotvel_staggered!`'s `lerp(H)` — and every one of them is evaluated through the same
# mask those kernels use, so an ice-free neighbour contributes exactly the zero it
# contributes to the residual itself.
#
# Row sum for the u-equation at `acx(i, j)`, writing P = ηH at `aa` and Q = hlerp(η)·lerp(H)
# at `ab` (derivation: expand ∂x(N_xx) + ∂y(N_xy) - βu and sum |coefficients| over both the
# u and the v unknowns, since the two equations are coupled):
#
#   ∂x(N_xx) : u-terms 8(P₋+P₊)/dx²   v-terms 4(P₋+P₊)/(dx·dy)
#   ∂y(N_xy) : u-terms 2(Q₋+Q₊)/dy²   v-terms 2(Q₋+Q₊)/(dx·dy)
#   basal    : β_face
#
# all divided by ρ·H_face. Sanity check: uniform η, H, no drag, dx = dy gives
# Λ = 32η/(ρdx²) hence Δτ = ρdx²/(16η) — which is Sandip's Eq. 7 at muB = 0, ndim2 = 4.1
# (ρdx²/16.4η). The two agree where they should; they differ only by the terms Eq. 7 omits.

@inline _etaH_aa(η, H, mask, i, j, k) =
    node_active(mask, NODE_AA, i, j) ? η[i, j, k] * H[i, j, k] : zero(eltype(η))

# `node_fully_active`, matching `_membrane_stress_staggered!`: `hlerp` of a zero-viscosity
# neighbour is `NaN`, not zero, so the strict rule is what makes this evaluable at all.
@inline _etaH_ab(η, H, mask, grid, i, j, k) =
    node_fully_active(mask, NODE_AB, i, j) ?
    hlerp(η, NODE_AB, grid, i, j, k) * lerp(H, NODE_AB, grid, i, j, k) : zero(eltype(η))

# Drag enters the bound only if the drag term actually depends on the velocity being
# iterated. Under `NoFrictionUpdate` the basal stress is a prescribed constant field, i.e. a
# forcing like the driving stress, contributing nothing to the operator's spectrum.
@inline _drag_in_spectrum(::ActiveFrictionUpdate) = true
@inline _drag_in_spectrum(::NoFrictionUpdate) = false

@kernel inbounds = true function _pseudo_dt_gershgorin!(dtau_x, dtau_y, η, H, β, ρ, scale,
                                                        drag, mask, dx, dy, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, k = I
    Z = zero(eltype(dtau_x))

    if node_active(mask, NODE_ACX, i, j)
        Hx = lerp(H, NODE_ACX, grid, I...)
        Ps = _etaH_aa(η, H, mask, i - 1, j, k) + _etaH_aa(η, H, mask, i, j, k)
        Qs = _etaH_ab(η, H, mask, grid, i, j, k) + _etaH_ab(η, H, mask, grid, i, j + 1, k)
        drag_x = drag ? lerp(β, NODE_ACX, grid, I...) : Z
        Λ = (8Ps / dx^2 + 4Ps / (dx * dy) + 2Qs / dy^2 + 2Qs / (dx * dy) + drag_x) / (ρ * Hx)
        dtau_x[I...] = (Hx > Z && Λ > Z) ? scale / Λ : Z
    else
        dtau_x[I...] = Z
    end

    if node_active(mask, NODE_ACY, i, j)
        Hy = lerp(H, NODE_ACY, grid, I...)
        Ps = _etaH_aa(η, H, mask, i, j - 1, k) + _etaH_aa(η, H, mask, i, j, k)
        Qs = _etaH_ab(η, H, mask, grid, i, j, k) + _etaH_ab(η, H, mask, grid, i + 1, j, k)
        drag_y = drag ? lerp(β, NODE_ACY, grid, I...) : Z
        Λ = (8Ps / dy^2 + 4Ps / (dx * dy) + 2Qs / dx^2 + 2Qs / (dx * dy) + drag_y) / (ρ * Hy)
        dtau_y[I...] = (Hy > Z && Λ > Z) ? scale / Λ : Z
    else
        dtau_y[I...] = Z
    end
end

"""
$(TYPEDSIGNATURES)

Fill `solver.dtau_x`/`dtau_y` for one solve, dispatching on `solver.pseudo_timestep`:
[`ViscosityPseudoTimeStep`](@ref) (the default) delegates to the Sandip Eq. 7 method above,
[`GershgorinPseudoTimeStep`](@ref) bounds the spectral radius of the residual operator the
solver actually iterates — including the basal-drag term and the ice mask (see the source
note above for the row sums, and [`GershgorinPseudoTimeStep`](@ref) for why they matter).

Called once by [`pseudo_transient!`](@ref) before the PT loop; `mask` and
`solver.friction_update` must be the same ones the loop itself will use, or the bound
describes a different operator than the one being iterated.
"""
function pseudo_dt!(solver::PseudoTransientSolver, mech::MechanicState, c::Constants,
                    rt::Runtime, mask::AbstractIceMask = NoMask())
    dx = Δx(rt.grid2d, Center(), 1, 1, 1)
    dy = Δy(rt.grid2d, Center(), 1, 1, 1)
    return pseudo_dt!(solver, solver.pseudo_timestep, mech, c, rt, mask, dx, dy)
end

function pseudo_dt!(solver::PseudoTransientSolver, pt::ViscosityPseudoTimeStep,
                    mech::MechanicState, c::Constants, rt::Runtime,
                    ::AbstractIceMask, dx, dy)
    pseudo_dt!(solver.dtau_x, solver.dtau_y, c.density_ice, dx, dy,
               mech.material.viscosity_depthaveraged, pt.muB, pt.ndim,
               solver.dtau_scaling, rt)
    return nothing
end

function pseudo_dt!(solver::PseudoTransientSolver, pt::GershgorinPseudoTimeStep,
                    mech::MechanicState, c::Constants, rt::Runtime,
                    mask::AbstractIceMask, dx, dy)
    gershgorin_dt!(solver, mech, c, rt, mask, 2 * pt.cfl * solver.dtau_scaling)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Fill `solver.dtau_x`/`dtau_y` with `scale / Λ`, `Λ` being the Gershgorin absolute row sum of
the mass-scaled residual operator (see the source note above). The prefactor is a parameter
because the two schemes reading this bound want different ones: the first-order iteration
takes `scale = 2·cfl` (`Δτ ≤ 2/λ_max`), [`AutotunedDynamicRelaxation`](@ref) takes
`scale = Δτ²` (`Δτ` enters the damped update twice, on the residual and on the accumulator).

`Λ = scale / dtau` is recoverable from the output, which is how the autotuner gets its
`M`-inner product without a third field.
"""
function gershgorin_dt!(solver::PseudoTransientSolver, mech::MechanicState, c::Constants,
                        rt::Runtime, mask::AbstractIceMask, scale)
    T = eltype(solver.dtau_x)
    dx = Δx(rt.grid2d, Center(), 1, 1, 1)
    dy = Δy(rt.grid2d, Center(), 1, 1, 1)
    rt.launch2d(rt.arch, rt.grid2d,
                _pseudo_dt_gershgorin! =>
                    (solver.dtau_x, solver.dtau_y, mech.material.viscosity_depthaveraged,
                     mech.topography.thickness, mech.friction.beta_eff,
                     convert(T, c.density_ice), convert(T, scale),
                     _drag_in_spectrum(solver.friction_update), mask,
                     convert(T, dx), convert(T, dy), rt.grid2d))
    return nothing
end

###############################################################
# Blatter-Pattyn Gershgorin pseudo-time step
###############################################################
#
# Row sums for BP's residual (`roadmaps/blatter-pattyn.md`, §2.1), per unit *volume* rather
# than per unit area: `P`/`Q` above are `ηH` at `aa`/`ab`; here they are bare `µ` (no `H`,
# matching `_membrane_stress_staggered_bp!`'s `σxx`/`σxy`), and there is a third,
# vertical-shear contribution `Λ_vert` the 2D bound has no counterpart for at all:
#
#   ∂x(σxx) : u-terms 8(P₋+P₊)/dx²   v-terms 4(P₋+P₊)/(dx·dy)
#   ∂y(σxy) : u-terms 2(Q₋+Q₊)/dy²   v-terms 2(Q₋+Q₊)/(dx·dy)
#   ∂z(σxz) : u-terms 2(R₋/δ₋ + R₊/δ₊)/Δ_k
#
# all divided by `ρ̃` (not `ρ̃H`) — `Λ_horiz` therefore carries **no** `H` at all (unlike the
# 2D bound's `/(ρHx)`), while `Λ_vert` alone carries the `1/H²` that falls out of `δ±`/`Δ_k`
# both being `H·Δζ` quantities (`R`/`µ` is already staggered in physical space via `hlerp`,
# but the two sigma-to-z conversions — one inside `σxz` itself via `x_dz`, one in the outer
# flux-difference — each contribute a `1/H`). At `k = 1` the bottom flux is the bed BC: its
# row-sum contribution is `β_face/(Δζ₁·Hx)`, **not** doubled and **not** `1/Hx²` — unlike an
# interior `R/δ` term, `β u_b` has no neighbouring layer to also appear as an off-diagonal
# for, and it carries only the single `1/H` `dotvel!`'s outer `_dz_over_H` division
# contributes (`β` itself, unlike `µ`, is not a per-length quantity). At `k = nz` the top
# flux is dropped outright (stress-free surface): no substitute term, matching `dotvel!`'s
# `top = k < nz ? sxz[k+1] : 0`. Sanity check to assert in tests: uniform `µ`, no drag,
# `dx = dy` and a single layer must reproduce the existing 2D bound divided by `H` — but
# `_check_momentum_grid` rejects `nz == 1` for BP outright, so this is only ever exercised at
# `nz ≥ 2`, where both `k = 1` and `k = nz` always have a genuine interior neighbour.

@inline _mu_aa(μ, mask, i, j, k) =
    node_active(mask, NODE_AA, i, j) ? μ[i, j, k] : zero(eltype(μ))

@inline _mu_ab(μ, mask, grid, i, j, k) =
    node_fully_active(mask, NODE_AB, i, j) ?
    hlerp(μ, NODE_AB, grid, i, j, k) : zero(eltype(μ))

# `µ` at the `acx_ac`/`acy_ac` interfaces `σxz`/`σyz` live on. **Not `node_fully_active`**,
# unlike `_mu_ab` above, and the asymmetry is forced rather than a preference:
# `_mask_cells` reads only the horizontal part of a node class, so `NODE_ACX_AC` and
# `NODE_ACX` resolve to the *same* cell pair `((i-1,j), (i,j))`. Gating `σxz` on the strict
# rule while the unknown itself is created under the loose one
# (`node_active(mask, NODE_ACX, ...)`, in `_dotvel_staggered_bp!`) therefore zeroes the entire
# vertical operator on every face of the margin ring — ~1 % of the solved faces on 8 km AIS.
# That is fatal for BP specifically: basal drag enters as the `k = 1` interface flux and
# nothing else couples the column to the bed, so layers `k ≥ 2` above such a face are left
# with no restoring force at all and drift linearly in pseudo-time
# (`roadmaps/blatter-pattyn.md`, Phase 1 "margin `σxz`"). `NODE_AB` is a genuine four-cell
# node, so the strict rule there stays right and stays put.
#
# One-sided instead of zero: interpolate down the ice-covered column alone (`NODE_AA_AC`,
# harmonic in `z` only), which is what `hlerp` onto `NODE_ACX_AC` already reduces to when the
# two columns carry equal `µ`. This never inverts an ice-free cell's `µ`, so it keeps the
# NaN-avoidance the strict check was reached for in the first place, and it is bit-for-bit
# the previous expression wherever both cells are active — i.e. everywhere but the margin.
@inline function _mu_acxz(μ, mask, grid, i, j, k)
    west = node_active(mask, NODE_AA, i - 1, j)
    east = node_active(mask, NODE_AA, i, j)
    return (west & east) ? hlerp(μ, NODE_ACX_AC, grid, i, j, k) :
           west ? hlerp(μ, NODE_AA_AC, grid, i - 1, j, k) :
           east ? hlerp(μ, NODE_AA_AC, grid, i, j, k) : zero(eltype(μ))
end

@inline function _mu_acyz(μ, mask, grid, i, j, k)
    south = node_active(mask, NODE_AA, i, j - 1)
    north = node_active(mask, NODE_AA, i, j)
    return (south & north) ? hlerp(μ, NODE_ACY_AC, grid, i, j, k) :
           south ? hlerp(μ, NODE_AA_AC, grid, i, j - 1, k) :
           north ? hlerp(μ, NODE_AA_AC, grid, i, j, k) : zero(eltype(μ))
end

# The three coefficient groups of §2.1, factored out because Phase 2's tridiagonal assembly
# (`_line_relax_x!`/`_line_relax_y!` below) must build the *same* operator this bound bounds —
# "one source of truth for the vertical operator, so the two cannot drift apart"
# (`roadmaps/blatter-pattyn.md`, Phase 2). All three are per unit volume and *before* the
# `1/ρ̃` mass scaling, which both callers apply themselves.

# `∂x(σxx) + ∂y(σxy)`, i.e. the membrane rows — the part that stays explicit under
# `ImplicitVertical` and therefore the part that alone sets `Δτ` there.
@inline function _lambda_horiz_x(μ, mask, grid, dx, dy, i, j, k)
    Ps = _mu_aa(μ, mask, i - 1, j, k) + _mu_aa(μ, mask, i, j, k)
    Qs = _mu_ab(μ, mask, grid, i, j, k) + _mu_ab(μ, mask, grid, i, j + 1, k)
    return 8Ps / dx^2 + 4Ps / (dx * dy) + 2Qs / dy^2 + 2Qs / (dx * dy)
end

@inline function _lambda_horiz_y(μ, mask, grid, dx, dy, i, j, k)
    Ps = _mu_aa(μ, mask, i, j - 1, k) + _mu_aa(μ, mask, i, j, k)
    Qs = _mu_ab(μ, mask, grid, i, j, k) + _mu_ab(μ, mask, grid, i + 1, j, k)
    return 8Ps / dy^2 + 4Ps / (dx * dy) + 2Qs / dx^2 + 2Qs / (dx * dy)
end

# One interface of `∂z(σxz)`: the magnitude of the `u[k∓1]` coupling of cell `k` across
# interface `kf` (`kf = k` for the bed-ward neighbour, `kf = k+1` for the surface-ward one),
# `R/(δ Δ_k H²)`. The Gershgorin row sum takes `2×` this (the coupling appears once on the
# diagonal and once off it); the tridiagonal takes it once each, on the two sides separately.
@inline _vshear_x(μ, mask, grid, Hx, i, j, k, kf) =
    _mu_acxz(μ, mask, grid, i, j, kf) / Δz(grid, Vertex(), i, j, kf) /
    Δz(grid, Center(), i, j, k) / Hx^2

@inline _vshear_y(μ, mask, grid, Hy, i, j, k, kf) =
    _mu_acyz(μ, mask, grid, i, j, kf) / Δz(grid, Vertex(), i, j, kf) /
    Δz(grid, Center(), i, j, k) / Hy^2

# The bed interface at `k = 1`, where the flux *is* `τ_b = β u_b` (§2.3): `β/(Δζ₁ H)`, with a
# single `1/H` and no factor 2 — unlike an interior `R/δ` term, `β u_b` has no neighbouring
# layer to also appear as an off-diagonal for, and `β` is not a per-length quantity.
@inline _vdrag_x(β, grid, grid2d, Hx, i, j) =
    lerp(β, NODE_ACX, grid2d, i, j, 1) / Δz(grid, Center(), i, j, 1) / Hx

@inline _vdrag_y(β, grid, grid2d, Hy, i, j) =
    lerp(β, NODE_ACY, grid2d, i, j, 1) / Δz(grid, Center(), i, j, 1) / Hy

# Whether the vertical rows belong in the *explicit* stability bound at all. Under
# `ImplicitVertical` they are inverted exactly rather than stepped over, so `Δτ` is bounded by
# `Λ_horiz` alone — the whole point of Phase 2, and the one line of the bound that changes.
@inline _bound_vertical(::ExplicitVertical, Λ_vert) = Λ_vert
@inline _bound_vertical(::ImplicitVertical, Λ_vert) = zero(Λ_vert)

@kernel inbounds = true function _pseudo_dt_gershgorin_bp!(dtau_x, dtau_y, μ, H, β, ρ, scale,
                                                            drag, nz, mask, dx, dy, grid,
                                                            grid2d, dtau_cap, vertical, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, k = I
    Z = zero(eltype(dtau_x))

    if node_active(mask, NODE_ACX, i, j)
        Hx = lerp(H, NODE_ACX, grid2d, i, j, 1)
        Λ_horiz = _lambda_horiz_x(μ, mask, grid, dx, dy, i, j, k)

        bottom = if k > 1
            2 * _vshear_x(μ, mask, grid, Hx, i, j, k, k)
        else
            drag ? _vdrag_x(β, grid, grid2d, Hx, i, j) : Z
        end
        top = k < nz ? 2 * _vshear_x(μ, mask, grid, Hx, i, j, k, k + 1) : Z

        Λ = Hx > Z ? (Λ_horiz + _bound_vertical(vertical, bottom + top)) / ρ : Z
        # `min(·, dtau_cap)`, not a mask-style ternary: a row whose vertical stiffness has
        # been zeroed by the strict `node_fully_active` check at a margin (bottom/top above)
        # still reports a mathematically valid — but locally much weaker — Λ than an interior
        # row of the same column, since the vertical term is normally dominant (§2.1). Left
        # uncapped (`dtau_cap = Inf`, the default), that is honest; a finite cap is an opt-in
        # numerical safety net for exactly that case, complementing (not replacing)
        # `clamp_velocity_gradients!`, which bounds the *consequence* rather than the cause.
        dtau_x[I...] = (Hx > Z && Λ > Z) ? min(scale / Λ, dtau_cap) : Z
    else
        dtau_x[I...] = Z
    end

    if node_active(mask, NODE_ACY, i, j)
        Hy = lerp(H, NODE_ACY, grid2d, i, j, 1)
        Λ_horiz = _lambda_horiz_y(μ, mask, grid, dx, dy, i, j, k)

        bottom = if k > 1
            2 * _vshear_y(μ, mask, grid, Hy, i, j, k, k)
        else
            drag ? _vdrag_y(β, grid, grid2d, Hy, i, j) : Z
        end
        top = k < nz ? 2 * _vshear_y(μ, mask, grid, Hy, i, j, k, k + 1) : Z

        Λ = Hy > Z ? (Λ_horiz + _bound_vertical(vertical, bottom + top)) / ρ : Z
        dtau_y[I...] = (Hy > Z && Λ > Z) ? min(scale / Λ, dtau_cap) : Z
    else
        dtau_y[I...] = Z
    end
end

"""
$(TYPEDSIGNATURES)

[`MomentumBalance3D`](@ref) counterpart of [`gershgorin_dt!`](@ref): fills
`solver.dtau_x`/`dtau_y` (`ACX3`/`ACY3`, one `Δτ` per layer, not just per column) with
`min(scale / Λ, dtau_cap)`, `Λ` being the row sum of §2.1 above. Same `scale` convention as
the 2D method: `2·cfl` for the plain iteration, `Δτ²` for [`AutotunedDynamicRelaxation`](@ref).

Which rows enter `Λ` is `solver.vertical_treatment`'s decision: under
[`ExplicitVertical`](@ref) all of them, under [`ImplicitVertical`](@ref) the membrane rows
alone, since the vertical operator is then inverted exactly rather than stepped over and no
longer constrains the explicit step.

`dtau_cap` (default `Inf`, a no-op) is a numerical safety net, not part of the row-sum
derivation: see the source note above `_pseudo_dt_gershgorin_bp!`.
"""
function gershgorin_dt!(solver::PseudoTransientSolver, mech::MechanicState, c::Constants,
                        rt::Runtime, momentum::MomentumBalance3D, mask::AbstractIceMask,
                        scale; dtau_cap = Inf)
    T = eltype(solver.dtau_x)
    dx = Δx(rt.grid2d, Center(), 1, 1, 1)
    dy = Δy(rt.grid2d, Center(), 1, 1, 1)
    nz = size(rt.grid, Center())[3]
    rt.launch(rt.arch, rt.grid,
              _pseudo_dt_gershgorin_bp! =>
                  (solver.dtau_x, solver.dtau_y, mech.material.viscosity,
                   mech.topography.thickness, mech.friction.beta_eff,
                   convert(T, c.density_ice), convert(T, scale),
                   _drag_in_spectrum(solver.friction_update), nz, mask,
                   convert(T, dx), convert(T, dy), rt.grid, rt.grid2d, convert(T, dtau_cap),
                   solver.vertical_treatment))
    return nothing
end

"""
$(TYPEDSIGNATURES)

Fill `solver.dtau_x`/`dtau_y` for one Blatter-Pattyn solve. Only
[`GershgorinPseudoTimeStep`](@ref) is implemented for [`MomentumBalance3D`](@ref) — Sandip's
Eq. 7 (`ViscosityPseudoTimeStep`) bounds only the membrane part of a *depth-integrated*
operator and has no vertical-shear term to extend, so it is not a meaningful bound for BP at
all, and passing it raises a `MethodError` rather than silently reusing the wrong formula.

`dtau_cap` (default `Inf`) is forwarded to [`gershgorin_dt!`](@ref).
"""
function pseudo_dt!(solver::PseudoTransientSolver, mech::MechanicState, c::Constants,
                    rt::Runtime, momentum::MomentumBalance3D, mask::AbstractIceMask = NoMask();
                    dtau_cap = Inf)
    dx = Δx(rt.grid2d, Center(), 1, 1, 1)
    dy = Δy(rt.grid2d, Center(), 1, 1, 1)
    return pseudo_dt!(solver, solver.pseudo_timestep, mech, c, rt, momentum, mask, dx, dy;
                      dtau_cap)
end

function pseudo_dt!(solver::PseudoTransientSolver, pt::GershgorinPseudoTimeStep,
                    mech::MechanicState, c::Constants, rt::Runtime,
                    momentum::MomentumBalance3D, mask::AbstractIceMask, dx, dy;
                    dtau_cap = Inf)
    gershgorin_dt!(solver, mech, c, rt, momentum, mask, 2 * pt.cfl * solver.dtau_scaling;
                   dtau_cap)
    return nothing
end

"""
$(TYPEDSIGNATURES)

The SSA/DIVA pseudo-transient velocity rate. The membrane-stress divergence is formed by Chmy's `∂x`/`∂y` acting directly on the
already-staggered `sxx` (`aa`), `sxy` (`ab`), `syy` (`aa`) — no interpolation, by the same
node-algebra argument as [`drivingstress!`](@ref) and [`velocitygradients!`](@ref): `∂x`
of `sxx` and `∂y` of `sxy` both land on `acx`; `∂x` of `sxy` and `∂y` of `syy` both land on
`acy`. Only the ice thickness `H` needs staggering onto the face, by `lerp`.

Depth-integrated throughout, so it runs on `rt.grid2d`.

The rate is accumulated with damping (Sandip et al. 2024, Eq. 12–14):
`dvx = (1 - gamma) * dvx_old + r_x(u)`, where `dvx_old` is whatever `dvx` already holds on
entry. `gamma = 1` (the default) discards `dvx_old` entirely and reproduces the plain,
undamped rate every call — this method's behaviour before `gamma` was added.
`0 < gamma < 1` gives the rate memory, which is what
lets [`pseudo_transient!`](@ref) converge in fewer iterations (the whole point of Sandip's
"acceleration owing to damping", their §2.5).

`resid_x`, `resid_y` are mandatory keyword outputs (no default: aliasing `dvx`/`dvy` would
corrupt the damped combination above, so a genuinely separate buffer is required) that
always receive the raw, undamped rate `r(u)/(ρH)` — independent of `gamma`, and thus a
meaningful "how far from momentum balance" diagnostic even when the accumulator itself is
damped. [`pseudo_transient!`](@ref) reduces these into the `residual` it returns.
"""
function dotvel!(dvx, dvy, sxx, sxy, syy, base_x, base_y, driving_x, driving_y,
                 H, density_ice, rt::Runtime,
                 momentum::MomentumBalance2D,
                 mask::AbstractIceMask = NoMask(); gamma = 1, resid_x, resid_y)
    ρ = convert(eltype(dvx), density_ice)
    γ = convert(eltype(dvx), gamma)
    rt.launch2d(rt.arch, rt.grid2d,
              _dotvel_staggered! => (dvx, dvy, resid_x, resid_y, sxx, sxy, syy, base_x,
                                     base_y, driving_x, driving_y, H, ρ, γ, mask, rt.grid2d))
    return nothing
end

@kernel inbounds = true function _dotvel_staggered!(dvx, dvy, resid_x, resid_y, sxx, sxy,
                                                     syy, base_x, base_y, driving_x,
                                                     driving_y, H, ρ, γ, mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    Z = zero(eltype(dvx))
    one_m_γ = one(γ) - γ

    if node_active(mask, NODE_ACX, i, j)
        Hx = lerp(H, NODE_ACX, grid, I...)
        shear_x = ∂x(sxx, grid, I...) + ∂y(sxy, grid, I...)
        raw_x = Hx > zero(Hx) ?
                (shear_x - base_x[I...] - driving_x[I...]) / (ρ * Hx) : Z
        resid_x[I...] = raw_x
        dvx[I...] = one_m_γ * dvx[I...] + raw_x
    else
        resid_x[I...] = Z
        dvx[I...] = Z
    end

    if node_active(mask, NODE_ACY, i, j)
        Hy = lerp(H, NODE_ACY, grid, I...)
        shear_y = ∂x(sxy, grid, I...) + ∂y(syy, grid, I...)
        raw_y = Hy > zero(Hy) ?
                (shear_y - base_y[I...] - driving_y[I...]) / (ρ * Hy) : Z
        resid_y[I...] = raw_y
        dvy[I...] = one_m_γ * dvy[I...] + raw_y
    else
        resid_y[I...] = Z
        dvy[I...] = Z
    end
end

###############################################################
# Blatter-Pattyn pseudo-transient velocity rate
###############################################################
#
# Three structural differences from `_dotvel_staggered!` above (`roadmaps/blatter-pattyn.md`
# §1, §2.3), all of them changes to what is assembled, not a new mechanism:
#
#  - a third divergence term, `∂z(σxz)`, on top of the same `∂x(σxx) + ∂y(σxy)` the 2D
#    kernel already computes (`sxx`/`sxy`/`syy` are per-unit-volume here, but the operator
#    algebra staggering them onto `acx`/`acy` is identical);
#  - **no `base_x`/`base_y` term in the residual body.** Basal drag has left the momentum
#    balance and become the *bed interface flux* inside `∂z(σxz)`'s own assembly (`bot`
#    below) — folding the BC into this kernel rather than a separate pass keeps `stress.xz`
#    honestly meaning `σxz` everywhere a caller might read it for diagnostics, never a
#    substituted boundary value;
#  - mass scaling is `/ρ̃` (a volume), not `/(ρ̃H)` (an area) — only the vertical term's
#    ζ→z conversion still needs `H`, via the same `_dz_over_H` `_velocity_gradients!` uses.
#
# `∂z(σxz)` is assembled as `_dz_over_H((top - bot) / Δζ_k, Hx)`, mirroring `∂z_σ`'s own
# `δ(f) * iΔ(...)` exactly (`Δz(grid, Center(), i, j, k)` is the same `Δζ_k` denominator
# `_viscosity_integrals!` already reads off the sigma axis) but with the two interface
# values substituted at the column ends per §2.3:
#
#   top = k < nz ? σxz[k+1] : 0          (surface: σxz = 0, stress-free)
#   bot = k > 1  ? σxz[k]   : τ_b,x      (bed: the flux *is* the basal traction)
#
# `τ_b,x` is `stress.base_x`, filled once per iteration by the same `update_basalstress!`
# SSA/DIVA already call — BP does not recompute it, it only relocates where it enters.

@kernel inbounds = true function _dotvel_staggered_bp!(dvx, dvy, resid_x, resid_y, sxx, sxy,
                                                        sxz, syy, syz, base_x, base_y,
                                                        driving_x, driving_y, H, ρ, γ, nz,
                                                        mask, grid, grid2d, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, k = I
    Z = zero(eltype(dvx))
    one_m_γ = one(γ) - γ

    if node_active(mask, NODE_ACX, i, j)
        Hx = lerp(H, NODE_ACX, grid2d, i, j, 1)
        dζ = Δz(grid, Center(), i, j, k)
        top = k < nz ? sxz[i, j, k + 1] : Z
        bot = k > 1  ? sxz[i, j, k]     : base_x[i, j, 1]
        vert_x = _dz_over_H((top - bot) / dζ, Hx)
        shear_x = ∂x(sxx, grid, I...) + ∂y(sxy, grid, I...) + vert_x
        raw_x = Hx > zero(Hx) ? (shear_x - driving_x[i, j, 1]) / ρ : Z
        resid_x[I...] = raw_x
        dvx[I...] = one_m_γ * dvx[I...] + raw_x
    else
        resid_x[I...] = Z
        dvx[I...] = Z
    end

    if node_active(mask, NODE_ACY, i, j)
        Hy = lerp(H, NODE_ACY, grid2d, i, j, 1)
        dζ = Δz(grid, Center(), i, j, k)
        top = k < nz ? syz[i, j, k + 1] : Z
        bot = k > 1  ? syz[i, j, k]     : base_y[i, j, 1]
        vert_y = _dz_over_H((top - bot) / dζ, Hy)
        shear_y = ∂x(sxy, grid, I...) + ∂y(syy, grid, I...) + vert_y
        raw_y = Hy > zero(Hy) ? (shear_y - driving_y[i, j, 1]) / ρ : Z
        resid_y[I...] = raw_y
        dvy[I...] = one_m_γ * dvy[I...] + raw_y
    else
        resid_y[I...] = Z
        dvy[I...] = Z
    end
end

"""
$(TYPEDSIGNATURES)

The Blatter-Pattyn pseudo-transient velocity rate — the [`dotvel!`](@ref) counterpart of
[`dotvel!(..., ::MomentumBalance2D, ...)`](@ref) for [`MomentumBalance3D`](@ref) — Robinson
et al. (2022) Eq. (1), assembled per unit volume with the bed/surface flux boundary
conditions folded into the vertical term (see the source note above).

Runs on `rt.grid`, the column grid. `resid_x`/`resid_y`, `gamma` and the damped-accumulator
semantics are otherwise identical to the 2D method: mandatory keyword outputs receiving the
raw, undamped rate, `gamma = 1` reproducing the plain rate every call.
"""
function dotvel!(dvx, dvy, sxx, sxy, sxz, syy, syz, base_x, base_y, driving_x, driving_y,
                 H, density_ice, rt::Runtime,
                 momentum::MomentumBalance3D,
                 mask::AbstractIceMask = NoMask(); gamma = 1, resid_x, resid_y)
    ρ = convert(eltype(dvx), density_ice)
    γ = convert(eltype(dvx), gamma)
    nz = size(rt.grid, Center())[3]
    rt.launch(rt.arch, rt.grid,
              _dotvel_staggered_bp! => (dvx, dvy, resid_x, resid_y, sxx, sxy, sxz, syy, syz,
                                        base_x, base_y, driving_x, driving_y, H, ρ, γ, nz,
                                        mask, rt.grid, rt.grid2d))
    return nothing
end

###############################################################
# Convergence criteria
###############################################################

# `_pt_error` is what the PT loop compares against `abstol`. Both methods reduce fields the
# loop already maintains, so neither costs an extra kernel — only the reduction the
# `ncheck` cadence exists to amortize.

@inline _abs_diff(a, b) = abs(a - b)

_pt_error(::VelocityIncrement, solver, ux, uy, ref) =
    max(mapreduce(_abs_diff, max, asarray(ux), asarray(solver.velocity_x_old)),
        mapreduce(_abs_diff, max, asarray(uy), asarray(solver.velocity_y_old)))

_pt_error(::ScaledResidual, solver, ux, uy, ref) =
    max(maximum(abs, asarray(solver.residual_x)),
        maximum(abs, asarray(solver.residual_y))) / ref

@kernel inbounds = true function _driving_rate!(rate_x, rate_y, driving_x, driving_y, H, ρ,
                                                 mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    Z = zero(eltype(rate_x))
    Hx = lerp(H, NODE_ACX, grid, I...)
    Hy = lerp(H, NODE_ACY, grid, I...)
    rate_x[I...] = (node_active(mask, NODE_ACX, i, j) && Hx > Z) ?
                   driving_x[I...] / (ρ * Hx) : Z
    rate_y[I...] = (node_active(mask, NODE_ACY, i, j) && Hy > Z) ?
                   driving_y[I...] / (ρ * Hy) : Z
end

# The scale `_pt_error(::ScaledResidual, ...)` divides by: the velocity rate the driving
# stress alone would produce, i.e. `dotvel!`'s residual with the membrane and basal terms
# dropped. Computed once per solve, before the loop.
#
# `solver.residual_x`/`residual_y` are borrowed as scratch: `dotvel!` overwrites them
# unconditionally on the first iteration, so nothing that is read afterwards is lost, and
# the alternative — two more `acx`/`acy` fields on every solver ever constructed — would be
# paid for by every solve to serve one convergence criterion.
_convergence_scale(::VelocityIncrement, mech, c, solver, rt, mask) =
    one(eltype(asarray(mech.velocity.depthaverage_x)))

function _convergence_scale(::ScaledResidual, mech::MechanicState, c::Constants,
                            solver::PseudoTransientSolver, rt::Runtime,
                            mask::AbstractIceMask)
    T = eltype(solver.residual_x)
    rt.launch2d(rt.arch, rt.grid2d,
                _driving_rate! => (solver.residual_x, solver.residual_y,
                                   mech.stress.driving_x, mech.stress.driving_y,
                                   mech.topography.thickness, convert(T, c.density_ice),
                                   mask, rt.grid2d))
    scale = max(maximum(abs, asarray(solver.residual_x)),
                maximum(abs, asarray(solver.residual_y)))
    # A domain with no driving stress at all (a flat, ice-free or perfectly level state) is
    # already in balance; falling back to 1 makes `err` the raw residual rather than `Inf`.
    return scale > zero(scale) ? scale : one(scale)
end

###############################################################
# Autotuned dynamic relaxation (Duretz et al. 2026)
###############################################################

"""
$(TYPEDSIGNATURES)

The iteration parameters in force at one point of the PT loop, as an immutable
`NamedTuple`: `gamma`, `theta_v`, the `Δτ` prefactor `scale` (`dtau = scale/Λ`, see
[`gershgorin_dt!`](@ref)), the last `lambda_min`, and the two partial sums
[`_arm_tuning`](@ref)/[`_tune!`](@ref) carry between iterations.

The field set is the same for every [`AbstractPTTuning`](@ref), so [`pseudo_transient!`](@ref)
stays type stable whichever is selected. `scale`/`lambda_min` are `NaN` under
[`FixedTuning`](@ref), which reads neither.
"""
_tuning_state(::Type{T}, gamma, theta_v, scale) where {T} =
    (; gamma = T(gamma), theta_v = T(theta_v), scale = T(scale), lambda_min = T(NaN),
       rayleigh_ur = zero(T), rayleigh_uu = zero(T), armed = false)

"""
$(TYPEDSIGNATURES)

Fill `solver.dtau_x`/`dtau_y` for the first iteration and return the initial
[`_tuning_state`](@ref).

[`FixedTuning`](@ref) defers to [`pseudo_dt!`](@ref) and freezes its own `gamma`/`theta_v`.
[`AutotunedDynamicRelaxation`](@ref) starts in its warm-up: undamped (`γ = 1`) at the
damped-Jacobi step `Δτ = 1/λ_max`, i.e. `scale = 1`.
"""
function _tuning_init!(tu::FixedTuning, solver::PseudoTransientSolver, mech::MechanicState,
                       c::Constants, rt::Runtime, mask::AbstractIceMask)
    pseudo_dt!(solver, mech, c, rt, mask)
    return _tuning_state(eltype(solver.dtau_x), tu.gamma, tu.theta_v, NaN)
end

function _tuning_init!(::AutotunedDynamicRelaxation, solver::PseudoTransientSolver,
                       mech::MechanicState, c::Constants, rt::Runtime,
                       mask::AbstractIceMask)
    T = eltype(solver.dtau_x)
    scale = T(solver.dtau_scaling)
    gershgorin_dt!(solver, mech, c, rt, mask, scale)
    return _tuning_state(T, 1, 1, scale)
end

@inline _du_dot_r(u, u_old, r) = (u - u_old) * r

# Guarded: `dtau` is exactly zero off-mask, where `Δu` is zero too, and an unguarded 0/0
# would `NaN` the whole reduction.
@inline _du2_over_dtau(u, u_old, dtau) =
    dtau > zero(dtau) ? (u - u_old)^2 / dtau : zero(dtau)

"""
$(TYPEDSIGNATURES)

`Δuᵀ r̃` over both velocity components, with `Δu = ux - solver.velocity_x_old`: the
numerator of the Rayleigh quotient, sampled at the two iterates that
[`_arm_tuning`](@ref) and [`_tune!`](@ref) each see one of.
"""
_sum_du_dot_r(solver, ux, uy) =
    mapreduce(_du_dot_r, +, asarray(ux), asarray(solver.velocity_x_old),
              asarray(solver.residual_x)) +
    mapreduce(_du_dot_r, +, asarray(uy), asarray(solver.velocity_y_old),
              asarray(solver.residual_y))

"""
$(TYPEDSIGNATURES)

`Σ Δu²/Δτ`, the Rayleigh quotient's `M`-inner product `Δuᵀ M Δu` up to the factor `scale`:
`M = diag(Λ)` is never stored, but [`gershgorin_dt!`](@ref) writes `dtau = scale/Λ`, so
`Λ = scale/dtau`.
"""
_sum_du2_over_dtau(solver, ux, uy) =
    mapreduce(_du2_over_dtau, +, asarray(ux), asarray(solver.velocity_x_old),
              asarray(solver.dtau_x)) +
    mapreduce(_du2_over_dtau, +, asarray(uy), asarray(solver.velocity_y_old),
              asarray(solver.dtau_y))

"""
$(TYPEDSIGNATURES)

Sample `Δuᵀ r̃(u^{k-1})` and `Δuᵀ M Δu`, and arm [`_tune!`](@ref) to take the missing
`Δuᵀ r̃(u^k)` next iteration. Splitting the Rayleigh quotient this way is what lets it be
evaluated from iterates the loop already holds, carrying two scalars instead of a second
residual buffer on every solver ever constructed.

Call after `pseudo_vel!`, where `ux`/`uy` hold `u^k`, `solver.velocity_*_old` hold
`u^{k-1}` and `solver.residual_*` hold `r̃(u^{k-1})`. A no-op under [`FixedTuning`](@ref),
and except every `cadence` iterations.
"""
_arm_tuning(::FixedTuning, state, solver, ux, uy, iter) = state

function _arm_tuning(tu::AutotunedDynamicRelaxation, state, solver::PseudoTransientSolver,
                     ux, uy, iter::Int)
    iter % tu.cadence == 0 || return state
    return merge(state, (; rayleigh_ur = _sum_du_dot_r(solver, ux, uy),
                           rayleigh_uu = _sum_du2_over_dtau(solver, ux, uy),
                           armed = true))
end

"""
$(TYPEDSIGNATURES)

Close the Rayleigh quotient armed by [`_arm_tuning`](@ref), derive `Δτ` and `γ` from it
(closed form in [`AutotunedDynamicRelaxation`](@ref)) and refill
`solver.dtau_x`/`dtau_y` — the refill doubling as the `λ_max` re-estimation a solve with an
evolving viscosity needs, since it rebuilds the Gershgorin bound from `η` as it now stands.

Call after `pseudo_rate!` and *before* the `u → u_old` copy, the one point where `Δu` is
still the increment [`_arm_tuning`](@ref) measured against and `solver.residual_*` already
hold `r̃(u^k)`. A no-op under [`FixedTuning`](@ref), and on a degenerate quotient (`Δu = 0`
at a stationary iterate), where the previous parameters stand.
"""
_tune!(::FixedTuning, state, solver, mech, c, rt, mask, ux, uy) = state

function _tune!(tu::AutotunedDynamicRelaxation, state, solver::PseudoTransientSolver,
                mech::MechanicState, c::Constants, rt::Runtime, mask::AbstractIceMask,
                ux, uy)
    state.armed || return state
    T = typeof(state.gamma)
    numerator = abs(_sum_du_dot_r(solver, ux, uy) - state.rayleigh_ur)
    denominator = state.scale * state.rayleigh_uu
    (numerator > 0 && denominator > 0) || return merge(state, (; armed = false))

    # A Rayleigh quotient of `Â` cannot exceed `λ_max ≤ 1`; above that it is measuring the
    # roundoff a converged `Δu` degenerates into, not the spectrum.
    λ_min = min(T(numerator / denominator), one(T))
    # `convert`ed, not promoted: `tuning` and `pseudo_timestep` are independent keyword
    # arguments, and a Float64 `cfl` on a Float32 solver must not retype `state`.
    d = 2 * convert(T, tu.c_damp) * sqrt(λ_min)
    cc = convert(T, solver.pseudo_timestep.cfl)^2
    dtau = -cc * d + sqrt(cc^2 * d^2 + 4 * cc)
    scale = dtau^2 * convert(T, solver.dtau_scaling)
    gershgorin_dt!(solver, mech, c, rt, mask, scale)
    return merge(state, (; gamma = d * dtau, scale, lambda_min = λ_min, armed = false))
end

###############################################################
# Glen viscosity continuation
###############################################################
#
# The dispatch types live in `src/mechanics/solvers.jl` next to `PseudoTransientSolver`;
# `update_viscosity!` lives here because it is PT-loop behaviour, called from
# `pseudo_rate!` below.

"""
$(TYPEDSIGNATURES)

No-op: `material.viscosity_depthaveraged` is untouched, since
[`NoViscosityContinuation`](@ref) means "treat viscosity as a fixed input".
"""
update_viscosity!(mech::MechanicState, ::NoViscosityContinuation, rt::Runtime,
                  mask::AbstractIceMask = NoMask()) = nothing

"""
$(TYPEDSIGNATURES)

[`GlenViscosityContinuation`](@ref): overwrite `material.viscosity_depthaveraged` in place
from the current velocity iterate — the effective strain rate
([`effective_strainrate_ssa!`](@ref), which requires `velocity`'s gradients to already be
current, i.e. [`velocitygradients!`](@ref) must run first) and `material.rate_factor_depthaveraged`,
relaxed in log-space toward the field's own previous value (see
[`GlenViscosityContinuation`](@ref)'s docstring for the formula).
"""
function update_viscosity!(mech::MechanicState, vc::GlenViscosityContinuation, rt::Runtime,
                          mask::AbstractIceMask = NoMask())
    (; velocity, strainrate, material) = mech
    effective_strainrate_ssa!(strainrate, velocity, rt, mask)
    rt.launch2d(rt.arch, rt.grid2d,
              _glen_viscosity_continuation! =>
                  (material.viscosity_depthaveraged, material.rate_factor_depthaveraged,
                   strainrate.effective_depthaveraged, vc.n_glen, vc.strainrate_reg,
                   vc.theta_mu, mask))
    return nothing
end

# No `grid` argument: every field here is already at `aa`, so unlike this file's other
# kernels there is no `lerp`/`hlerp` staggering to do.
@kernel inbounds = true function _glen_viscosity_continuation!(μ, A, eff, n, ε̇0, θ, mask, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    if node_active(mask, NODE_AA, i, j)
        ε̇ = sqrt(eff[I...]^2 + ε̇0^2)
        μ_raw = inv(2 * A[I...]^(1 / n)) * ε̇^((1 - n) / n)
        μ[I...] = exp(θ * log(μ_raw) + (1 - θ) * log(μ[I...]))
    end
end

"""
$(TYPEDSIGNATURES)

[`DIVAViscosityContinuation`](@ref): the same Glen law as the depth-averaged method above,
evaluated **per layer** on the column grid. Writes `material.viscosity` (`µ(z)`) from
DIVA's effective strain rate ([`effective_strainrate_diva!`](@ref), Eq. 13) and the column
rate factor `material.rate_factor`, then derives `material.viscosity_depthaveraged` from it
by [`depthaverage!`](@ref).

`µ̄` is *derived*, never computed independently: the membrane stress and the vertical shear
must describe the same ice, so `µ̄` has to be the average of exactly the `µ(z)` the shear was
built from.

Reuses `_glen_viscosity_continuation!` unchanged — the kernel indexes every field at its own
`I`, so the identical code serves the `aa`/`grid2d` and `aa`/`grid` launches.
"""
function update_viscosity!(mech::MechanicState, vc::DIVAViscosityContinuation, rt::Runtime,
                           mask::AbstractIceMask = NoMask())
    (; material, strainrate) = mech
    effective_strainrate_diva!(mech, rt, mask)
    rt.launch(rt.arch, rt.grid,
              _glen_viscosity_continuation! =>
                  (material.viscosity, material.rate_factor, strainrate.effective,
                   vc.n_glen, vc.strainrate_reg, vc.theta_mu, mask))
    depthaverage!(material.viscosity_depthaveraged, material.viscosity, rt, mask)
    return nothing
end

"""
$(TYPEDSIGNATURES)

[`BPViscosityContinuation`](@ref): the same per-layer Glen law as
[`DIVAViscosityContinuation`](@ref)'s method above, but from BP's effective strain rate
([`effective_strainrate_bp!`](@ref), Eq. 3, built from the real `velocity.x_dz`/`y_dz`
rather than diagnosed from `τ_b`) and reusing `_glen_viscosity_continuation!` unchanged.

Does **not** derive `material.viscosity_depthaveraged`: nothing on the BP path reads it
(`roadmaps/blatter-pattyn.md`, §1.3), so — unlike the DIVA method — there is no
`depthaverage!` call here, and the field is left at whatever degenerate allocation the state
constructor gave it.
"""
function update_viscosity!(mech::MechanicState, vc::BPViscosityContinuation, rt::Runtime,
                           mask::AbstractIceMask = NoMask())
    (; material, strainrate) = mech
    effective_strainrate_bp!(mech, rt, mask)
    rt.launch(rt.arch, rt.grid,
              _glen_viscosity_continuation! =>
                  (material.viscosity, material.rate_factor, strainrate.effective,
                   vc.n_glen, vc.strainrate_reg, vc.theta_mu, mask))
    return nothing
end

###############################################################
# DIVA depth-integrated-viscosity chain
###############################################################
#
# The full chain, in dependency order:
#
#   µ(z)  ←  ε̇_e (Eq. 13, u_z from Eq. 21 using the *previous* τ_b and µ)
#   µ̄     ←  depthaverage(µ(z))                                   [both in update_viscosity!]
#   F₁,F₂ ←  ∫(1/µ)((s-z)/H)^m dz   (Eq. 15)
#   β_eff ←  1/(1/β + F₂)           (Eqs. 19–20)
#
# `τ_b` is deliberately *not* recomputed here: it is read as whatever the last
# `update_basalstress!` left, which is exactly the paper's "obtained from the previous
# iteration". On a cold start it is zero, so `u_z = 0` and the first `ε̇_e` is the SSA
# invariant — a sane starting point rather than a singular one.

"""
$(TYPEDSIGNATURES)

Evaluate DIVA's depth-integrated-viscosity chain once, in dependency order: `µ(z)` and `µ̄`
via [`update_viscosity!`](@ref), then `F₁`/`F₂` via [`viscosity_integrals!`](@ref), then
`β_eff` via [`beta_eff_diva!`](@ref).

This is the function a caller must invoke before [`pseudo_transient!`](@ref) under the
default [`NoDIVUpdate`](@ref) — the solver will not do it for you (`roadmaps/chmy.md`,
Phase 3, decision 7). Under [`PeriodicDIVUpdate`](@ref) the loop also calls it every
`n_update` iterations.

!!! warning "Skipping this leaves `β_eff = 0`, i.e. frictionless sliding"
    `friction.beta_eff` is zero at allocation, and a zero friction coefficient is a
    perfectly well-formed (if physically absurd) input that raises no error — the solve will
    simply run away. There is deliberately no runtime guard; this docstring is the contract,
    exactly as with the halo-filling requirement on [`pseudo_transient!`](@ref).

Requires `mech.material.viscosity` to hold a usable previous iterate (the chain's `ε̇_e`
divides by it) and the depth-averaged velocity gradients to be current.
"""
function diva_update!(mech::MechanicState, solver::PseudoTransientSolver, rt::Runtime,
                      mask::AbstractIceMask = NoMask())
    (; material) = mech
    update_viscosity!(mech, solver.viscosity_continuation, rt, mask)
    viscosity_integrals!(material.viscosity_integral_1, material.viscosity_integral_2,
                         mech, rt, mask)
    beta_eff_diva!(mech, rt, mask)
    return nothing
end

# "How often", per decision 6 — kept strictly separate from "how" (the continuation type).
# SSA has no DIV chain, so its viscosity continuation keeps its every-iteration cadence
# inside `pseudo_rate!` untouched; only the DIVA path consults `solver.div_update`.
_div_refresh_due(::NoDIVUpdate, iter) = false
_div_refresh_due(d::PeriodicDIVUpdate, iter) = iter % d.n_update == 0

# Who drives the viscosity from inside `pseudo_rate!`, i.e. every PT iteration.
#
# SSA: its continuation (`GlenViscosityContinuation`) writes `viscosity_depthaveraged`
# directly from Eq. 12, and doing so every iteration *is* the continuation — that is what
# lets the nonlinear solve relax rather than jump. Unchanged.
#
# DIVA: nothing here. Its continuation writes `viscosity` (`µ(z)`) and only `diva_update!`
# may drive it, at the cadence `solver.div_update` names. Calling it here as well would
# refresh `µ(z)` every iteration regardless of that strategy — the exact contradiction the
# "how" / "how often" split exists to prevent. The two continuations write *different*
# fields (`AA2` vs `AA3`), which is why SSA's cadence needs no gate at all.
_iterate_viscosity!(mech::MechanicState, ::SSAMomentumBalance,
                    solver::PseudoTransientSolver, rt::Runtime, mask::AbstractIceMask) =
    update_viscosity!(mech, solver.viscosity_continuation, rt, mask)

_iterate_viscosity!(::MechanicState, ::DIVAMomentumBalance, ::PseudoTransientSolver,
                    ::Runtime, ::AbstractIceMask) = nothing

# BP: like SSA, not DIVA. `BPViscosityContinuation` writes `viscosity` (`µ(z)`) directly
# from the real velocity gradients every iteration — there is no `τ_b`/`F₂` fixed point to
# stagger it against (`roadmaps/blatter-pattyn.md`, §2.2: "no F₁/F₂ chain and no β_eff — BP
# has no depth-integrated closure to build"), so `solver.div_update` plays no role here.
_iterate_viscosity!(mech::MechanicState, ::BlatterPattynMomentumBalance,
                    solver::PseudoTransientSolver, rt::Runtime, mask::AbstractIceMask) =
    update_viscosity!(mech, solver.viscosity_continuation, rt, mask)

###############################################################
# Basal friction update
###############################################################
#
# Split between `src/mechanics/solvers.jl` and here for the same reason as the
# viscosity-continuation trio above.

"""
$(TYPEDSIGNATURES)

[`ActiveFrictionUpdate`](@ref): the ordinary basal friction law. Overwrites
`velocity.base_x`/`base_y` from the current (SSA-limit) depth-averaged velocity, then
`stress.base_x`/`base_y` from `friction.beta_eff * velocity.base_{x,y}` via
[`basalstress!`](@ref) — the behaviour every solver had before
[`AbstractFrictionUpdate`](@ref) existed.
"""
function update_basalstress!(mech::MechanicState, ::ActiveFrictionUpdate, rt::Runtime,
                             mask::AbstractIceMask = NoMask())
    (; velocity, stress, friction) = mech
    # The SSA limit `u_b = ū`, which is what both balances use today. DIVA's correction
    # `u_b = ū/(1 + βF₂)` (Robinson et al. 2022, Eq. 18) is Stage 2 — see
    # `roadmaps/chmy.md`, Phase 3. Reads `depthaverage_x`/`y`, the field the solver
    # actually iterates (decision 1); it used to read `velocity.x`/`y`, which held the same
    # numbers only because `nz == 1` collapsed `ACX3` onto `ACX2`.
    copyto!(asarray(velocity.base_x), asarray(velocity.depthaverage_x))
    copyto!(asarray(velocity.base_y), asarray(velocity.depthaverage_y))
    basalstress!(stress.base_x, stress.base_y, friction.beta_eff,
                velocity.base_x, velocity.base_y, rt, mask)
    return nothing
end

"""
$(TYPEDSIGNATURES)

No-op: `stress.base_x`/`base_y` are untouched, since [`NoFrictionUpdate`](@ref) means
"hold the basal stress fixed at whatever was written there before the solve" — bypassing
the friction law entirely.
"""
update_basalstress!(::MechanicState, ::NoFrictionUpdate, ::Runtime,
                    ::AbstractIceMask = NoMask()) = nothing

"""
$(TYPEDSIGNATURES)

Evaluate the PT velocity rate: velocity gradients, viscosity continuation
([`update_viscosity!`](@ref), a no-op unless `solver.viscosity_continuation` enables it —
must run after the gradients it reads and before the membrane stress that reads its
output, hence its place in this order), membrane stress, basal stress
([`update_basalstress!`](@ref), a no-op if `solver.friction_update` is
[`NoFrictionUpdate`](@ref) — from the SSA-limit basal velocity otherwise; a real DIVA
`F₂` integral is Phase 3 future work), then the PT
velocity rate. Writes `solver.velocity_x_dt`/`velocity_y_dt`, damped by `gamma`, and
`solver.residual_x`/`residual_y`, the undamped rate (see [`dotvel!`](@ref)).

`gamma` is a keyword because [`AutotunedDynamicRelaxation`](@ref) recomputes the damping
*during* the solve, so the value in force at a given iteration is not a property of the
solver; it defaults to `1` (undamped), matching [`dotvel!`](@ref).
"""
# Per-balance grid requirements. SSA is depth-independent by construction: it never reads a
# column field, so any `nz` is fine and a column grid simply costs nothing. DIVA is the
# opposite — it is *defined* by its vertical structure, and on `nz == 1` the single
# quadrature point sits at `σ = ½`, which makes `F₂ = H/(4µ)` against the true `H/(3µ)`
# (see `viscosity_integrals!`). That is 25% low while looking entirely plausible, so it is
# an error rather than a warning (`roadmaps/chmy.md`, Phase 3, decision 4).
_check_momentum_grid(::SSAMomentumBalance, ::Runtime) = nothing

function _check_momentum_grid(::DIVAMomentumBalance, rt::Runtime)
    rt.grid === rt.grid2d && throw(ArgumentError(
        "DIVAMomentumBalance requires a column StaggeredGrid (nz > 1). On a " *
        "depth-averaged grid (nz == 1, rt.grid === rt.grid2d) the viscosity integral F₂ " *
        "is 25% low — a plausible-looking wrong answer, not an approximation. Build the " *
        "grid with a `layering` argument, e.g. " *
        "`StaggeredGrid(T, lx, ly, dx, dy, CorrectedVerticalLayering(T, " *
        "QuadraticSigmaTransform(T, nz)))`, or use SSAMomentumBalance()."))
    return nothing
end

# BP is stricter still than DIVA: on `nz == 1` it does not merely lose 25% accuracy in a
# derived quantity, it silently *degenerates to SSA* — the vertical-shear divergence
# `∂z(µ u_z)` a single layer cannot resolve simply vanishes, and the remaining membrane-only
# balance is a plausible-looking wrong answer with the wrong (3D, unaveraged) viscosity to
# boot (`roadmaps/blatter-pattyn.md`, Phase 1).
function _check_momentum_grid(::BlatterPattynMomentumBalance, rt::Runtime)
    rt.grid === rt.grid2d && throw(ArgumentError(
        "BlatterPattynMomentumBalance requires a column StaggeredGrid (nz > 1). On a " *
        "depth-averaged grid (nz == 1, rt.grid === rt.grid2d) the vertical-shear term " *
        "∂z(µ u_z) cannot be resolved and the solve silently degenerates to SSA with the " *
        "wrong (3D, unaveraged) viscosity — a plausible-looking wrong answer, not an " *
        "approximation. Build the grid with a `layering` argument, e.g. " *
        "`StaggeredGrid(T, lx, ly, dx, dy, CorrectedVerticalLayering(T, " *
        "QuadraticSigmaTransform(T, nz)))`, or use SSAMomentumBalance()."))
    return nothing
end

function pseudo_rate!(mech::MechanicState, c::Constants, rt::Runtime,
                      momentum::MomentumBalance2D,
                      solver::PseudoTransientSolver, mask::AbstractIceMask = NoMask();
                      gamma = 1)
    (; velocity, material, topography, stress) = mech

    depthaverage_velocitygradients!(velocity, rt, mask)
    _iterate_viscosity!(mech, momentum, solver, rt, mask)
    membranestress!(stress, velocity, material, topography, momentum, rt, mask)

    update_basalstress!(mech, solver.friction_update, rt, mask)

    dotvel!(solver.velocity_x_dt, solver.velocity_y_dt,
           stress.membrane_xx, stress.membrane_xy, stress.membrane_yy,
           stress.base_x, stress.base_y, stress.driving_x, stress.driving_y,
           topography.thickness, c.density_ice, rt, momentum, mask;
           gamma, resid_x = solver.residual_x, resid_y = solver.residual_y)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Iterate `mech.velocity.x`/`y` in pseudo-time until the SSA/DIVA momentum-balance residual
vanishes. Returns a named tuple
`(; iterations, error, converged, residual, damping, lambda_min)`. `error`/`converged` are
decided by `solver.convergence`; `residual` is the max-norm raw rate `‖r(u)/(ρH)‖` (`solver.residual_x`/`residual_y`, see
[`dotvel!`](@ref)) at the final iterate — always reported, whether or not it also gates the
loop; `damping` is the `γ` in force at the end of the solve and `lambda_min` the last
spectral estimate behind it (`lambda_min` is `NaN` unless `solver.tuning` is
[`AutotunedDynamicRelaxation`](@ref)).

Three solver-held strategies steer the iteration:

 - `solver.pseudo_timestep` ([`AbstractPseudoTimeStep`](@ref)) — how `Δτ` is chosen.
   [`GershgorinPseudoTimeStep`](@ref) (default) bounds the residual operator's spectral
   radius, basal drag included; [`ViscosityPseudoTimeStep`](@ref) is Sandip's `Δτ ∝ 1/η`,
   which **diverges on real geometry with a friction law being solved for** because it
   omits `β/(ρH)`.
 - `solver.convergence` ([`AbstractPTConvergence`](@ref)) — what `abstol` measures.
   [`VelocityIncrement`](@ref) (default) is the exact max-norm increment `|u_new - u_old|`
   (`ux`/`uy` still hold the pre-update iterate in `ux_old`/`uy_old` at that point) rather
   than the algebraic reconstruction `theta_v * dtau * dv` an earlier version used — that
   reconstruction relied on `dtau` being the same scalar everywhere, which a local `Δτ`
   field is not. [`ScaledResidual`](@ref) is the driving-stress-normalized momentum
   residual, and is what a domain mixing drag regimes (grounded ice + ice shelves) needs:
   the increment criterion cannot distinguish a converged shelf from a stalled one.
 - `solver.tuning` ([`AbstractPTTuning`](@ref)) — where `Δτ` and the damping come from.
   [`FixedTuning`](@ref) (default) carries them as hand-set numbers;
   [`AutotunedDynamicRelaxation`](@ref) derives both from a Gershgorin `λ_max` and a
   Rayleigh-quotient `λ_min`, re-derived every `cadence` iterations.

!!! warning "`material.viscosity_depthaveraged`, `topography.thickness` and `friction.beta_eff` need their halo filled by the caller"
    Unlike the velocity, these are fixed inputs for the whole solve, so this function does
    not refresh their halo itself (that is a Phase 4, `bc!`-policy decision belonging to
    whatever step last wrote them, not to the solver reading them). Leaving a halo at its
    allocation default (`setdata!` only touches the interior, by design — see its
    docstring) is not merely inexact at the domain edge, it can be `NaN`: the membrane
    stress's `hlerp(η)` reads the outer ghost ring of `viscosity_depthaveraged` at the two
    boundary `ab` corners, and `hlerp` of an unset (zero) neighbour is `NaN`, exactly like
    the already-documented ice-free case — except here there is no ice-free cell at all,
    only an unfilled halo. Found while writing this solver's own tests, which fill every
    such field with `fill_analytic!` for exactly this reason.

Requires a depth-averaged grid (`rt.grid2d === rt.grid`, i.e. `nz == 1`, checked): DIVA's
vertical-shear integral is Phase 3 future work, so today `mech.velocity.x`/`y` *are* the
depth-averaged velocity being solved for.

Two accelerations over the plain Sandip iteration (`roadmaps/PT-autotune.md` Phase 1):
the pseudo-time step is a *local* field (`solver.dtau_x`/`dtau_y`, via
[`pseudo_dt!`](@ref)) rather than one scalar throttled by the domain's single stiffest
cell, and the velocity rate is damped (via [`dotvel!`](@ref)) rather than recomputed from
scratch every iteration, which is what buys sub-quadratic iteration scaling.
"""
function pseudo_transient!(mech::MechanicState, c::Constants, solver::PseudoTransientSolver,
                           rt::Runtime,
                           momentum::MomentumBalance2D = SSAMomentumBalance(),
                           mask::AbstractIceMask = NoMask())
    _check_momentum_grid(momentum, rt)

    (; velocity) = mech
    (; abstol, maxiter, printout_every, ncheck, tuning) = solver

    ux, uy = velocity.depthaverage_x, velocity.depthaverage_y
    ux_old, uy_old = solver.velocity_x_old, solver.velocity_y_old
    dvx, dvy = solver.velocity_x_dt, solver.velocity_y_dt
    dtau_x, dtau_y = solver.dtau_x, solver.dtau_y
    resid_x, resid_y = solver.residual_x, solver.residual_y

    # Fields held fixed over the PT iteration.
    drivingstress!(mech, c, rt, mask)
    state = _tuning_init!(tuning, solver, mech, c, rt, mask)
    # After `drivingstress!` (it reads the driving stress) and before the loop (it borrows
    # `residual_x`/`residual_y` as scratch, which `dotvel!` overwrites on iteration 1).
    scale = _convergence_scale(solver.convergence, mech, c, solver, rt, mask)

    T = eltype(asarray(ux))
    err  = typemax(T)
    iter = 0
    while err > abstol && iter < maxiter
        iter += 1

        # DIVA's depth-integrated-viscosity chain, at the cadence `solver.div_update` names
        # (never, under the default `NoDIVUpdate`). Before `pseudo_rate!`, so the iteration
        # sees the refreshed `β_eff`/`µ̄` rather than applying them one iteration late.
        #
        # `pseudo_dt!` follows it because `β_eff` feeds the Gershgorin bound behind
        # `dtau_x`/`dtau_y`: a `β_eff` that grew under a stale bound makes that bound
        # optimistic, and the explicit iteration then diverges (decision 8).
        if momentum isa DIVAMomentumBalance && _div_refresh_due(solver.div_update, iter)
            diva_update!(mech, solver, rt, mask)
            pseudo_dt!(solver, mech, c, rt, mask)
        end

        pseudo_rate!(mech, c, rt, momentum, solver, mask; gamma = state.gamma)
        state = _tune!(tuning, state, solver, mech, c, rt, mask, ux, uy)

        # After `_tune!`, which needs the pre-copy `u_old` (see its docstring).
        copyto!(asarray(ux_old), asarray(ux))
        copyto!(asarray(uy_old), asarray(uy))

        pseudo_vel!(asarray(ux), asarray(ux_old), asarray(dvx), asarray(dtau_x),
                    state.theta_v)
        pseudo_vel!(asarray(uy), asarray(uy_old), asarray(dvy), asarray(dtau_y),
                    state.theta_v)

        bc!(rt.arch, rt.grid2d, ux => Neumann())
        bc!(rt.arch, rt.grid2d, uy => Neumann())

        state = _arm_tuning(tuning, state, solver, ux, uy, iter)

        if iter % ncheck == 0 || iter == maxiter
            err = _pt_error(solver.convergence, solver, ux, uy, scale)
        end
        if iter % printout_every == 0
            println("PT iteration $iter: err = $err, gamma = $(state.gamma)")
        end
    end

    # Diagnostic only (does not gate the loop above): the raw, undamped rate norm at the
    # final iterate, independent of `gamma` — see `dotvel!`'s `resid_x`/`resid_y` note.
    residual = max(maximum(abs, asarray(resid_x)), maximum(abs, asarray(resid_y)))

    return (; iterations = iter, error = err, converged = err <= abstol, residual,
              damping = state.gamma, lambda_min = state.lambda_min)
end

###############################################################
# Blatter-Pattyn pseudo-transient solve
###############################################################
#
# `roadmaps/blatter-pattyn.md`, Phase 1. Every piece below is a new method distinguished
# from its `MomentumBalance2D` counterpart either by dispatching on the disjoint
# `MomentumBalance3D` union or, where the existing function takes no `momentum` argument at
# all (`update_basalstress!`, `pseudo_dt!`, `gershgorin_dt!`, `_tuning_init!`, `_tune!`,
# `_convergence_scale`), by a new method one argument longer — so nothing here can ever be
# reached by a 2D call, and no 2D method body is touched.

# `_arm_tuning` is not repeated here: it already takes no `momentum`/grid argument at all
# (`solver.residual_x`/`velocity_x_old`/`dtau_x` alone, reduced via `asarray`), so the
# existing method serves BP unchanged (`roadmaps/blatter-pattyn.md`, §1.2).

@kernel inbounds = true function _basal_velocity_bp!(base_x, base_y, x, y, mask, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    base_x[i, j, 1] = node_active(mask, NODE_ACX, i, j) ? x[i, j, 1] : zero(eltype(base_x))
    base_y[i, j, 1] = node_active(mask, NODE_ACY, i, j) ? y[i, j, 1] : zero(eltype(base_y))
end

"""
$(TYPEDSIGNATURES)

[`MomentumBalance3D`](@ref) counterpart of
[`update_basalstress!(..., ::ActiveFrictionUpdate, ::Runtime)`](@ref): overwrites
`velocity.base_x`/`base_y` from the **actual** basal-layer velocity `velocity.x`/`y[:,:,1]`
— not the SSA-limit `depthaverage_x`/`y` the 2D method reads, which BP never populates (a
degenerate `1×1` allocation on this path, `roadmaps/blatter-pattyn.md`, §1.3) — then
`stress.base_x`/`base_y` from `friction.beta_eff * velocity.base_{x,y}` via
[`basalstress!`](@ref), exactly as the 2D method does downstream.
"""
function update_basalstress!(mech::MechanicState, ::ActiveFrictionUpdate, rt::Runtime,
                             ::MomentumBalance3D, mask::AbstractIceMask = NoMask())
    (; velocity, stress, friction) = mech
    rt.launch2d(rt.arch, rt.grid2d,
                _basal_velocity_bp! => (velocity.base_x, velocity.base_y, velocity.x,
                                        velocity.y, mask))
    basalstress!(stress.base_x, stress.base_y, friction.beta_eff,
                velocity.base_x, velocity.base_y, rt, mask)
    return nothing
end

update_basalstress!(::MechanicState, ::NoFrictionUpdate, ::Runtime, ::MomentumBalance3D,
                    ::AbstractIceMask = NoMask()) = nothing

"""
$(TYPEDSIGNATURES)

Evaluate the Blatter-Pattyn PT velocity rate: column velocity gradients
([`velocitygradients!`](@ref), then [`clamp_velocity_gradients!`](@ref) if `strainrate_cap`
is finite), viscosity continuation, membrane stress ([`membranestress!`](@ref)), basal
stress, then [`dotvel!`](@ref) — the same order as
[`pseudo_rate!(..., ::MomentumBalance2D, ...)`](@ref), minus the DIVA-only
`diva_update!` step BP has no counterpart for (`roadmaps/blatter-pattyn.md`, §2.2).

`strainrate_cap` (default `Inf`, a no-op) is the numerical safety net documented at
[`clamp_velocity_gradients!`](@ref).
"""
function pseudo_rate!(mech::MechanicState, c::Constants, rt::Runtime,
                      momentum::MomentumBalance3D,
                      solver::PseudoTransientSolver, mask::AbstractIceMask = NoMask();
                      gamma = 1, strainrate_cap = Inf)
    (; velocity, material, stress, topography) = mech

    velocitygradients!(velocity, topography.thickness, rt, mask)
    clamp_velocity_gradients!(velocity, strainrate_cap, rt)
    _iterate_viscosity!(mech, momentum, solver, rt, mask)
    membranestress!(stress, velocity, material, momentum, rt, mask)

    update_basalstress!(mech, solver.friction_update, rt, momentum, mask)

    dotvel!(solver.velocity_x_dt, solver.velocity_y_dt,
           stress.xx, stress.xy, stress.xz, stress.yy, stress.yz,
           stress.base_x, stress.base_y, stress.driving_x, stress.driving_y,
           topography.thickness, c.density_ice, rt, momentum, mask;
           gamma, resid_x = solver.residual_x, resid_y = solver.residual_y)
    return nothing
end

function _tuning_init!(tu::FixedTuning, solver::PseudoTransientSolver, mech::MechanicState,
                       c::Constants, rt::Runtime, momentum::MomentumBalance3D,
                       mask::AbstractIceMask; dtau_cap = Inf)
    pseudo_dt!(solver, mech, c, rt, momentum, mask; dtau_cap)
    return _tuning_state(eltype(solver.dtau_x), tu.gamma, tu.theta_v, NaN)
end

function _tuning_init!(::AutotunedDynamicRelaxation, solver::PseudoTransientSolver,
                       mech::MechanicState, c::Constants, rt::Runtime,
                       momentum::MomentumBalance3D, mask::AbstractIceMask; dtau_cap = Inf)
    T = eltype(solver.dtau_x)
    scale = T(solver.dtau_scaling)
    gershgorin_dt!(solver, mech, c, rt, momentum, mask, scale; dtau_cap)
    return _tuning_state(T, 1, 1, scale)
end

_tune!(::FixedTuning, state, solver, mech, c, rt, ::MomentumBalance3D, mask, ux, uy;
      dtau_cap = Inf) = state

function _tune!(tu::AutotunedDynamicRelaxation, state, solver::PseudoTransientSolver,
                mech::MechanicState, c::Constants, rt::Runtime, momentum::MomentumBalance3D,
                mask::AbstractIceMask, ux, uy; dtau_cap = Inf)
    state.armed || return state
    T = typeof(state.gamma)
    numerator = abs(_sum_du_dot_r(solver, ux, uy) - state.rayleigh_ur)
    denominator = state.scale * state.rayleigh_uu
    (numerator > 0 && denominator > 0) || return merge(state, (; armed = false))
    λ_min = min(T(numerator / denominator), one(T))
    d = 2 * convert(T, tu.c_damp) * sqrt(λ_min)
    cc = convert(T, solver.pseudo_timestep.cfl)^2
    dtau = -cc * d + sqrt(cc^2 * d^2 + 4 * cc)
    scale = dtau^2 * convert(T, solver.dtau_scaling)
    gershgorin_dt!(solver, mech, c, rt, momentum, mask, scale; dtau_cap)
    return merge(state, (; gamma = d * dtau, scale, lambda_min = λ_min, armed = false))
end

# `VelocityIncrement` needs no field to normalize against on either grid, so the 3D method
# just avoids the 2D one's (harmless but pointless) read of the degenerate `depthaverage_x`
# for its `eltype`, reading a field BP actually owns instead.
_convergence_scale(::VelocityIncrement, mech, c, solver::PseudoTransientSolver, rt,
                   ::MomentumBalance3D, mask) =
    one(eltype(asarray(solver.velocity_x_old)))

# `ScaledResidual`'s 2D implementation launches `_driving_rate!` on `rt.grid2d` into
# `solver.residual_x`/`residual_y` — for BP those are column-shaped (`ACX3`/`ACY3`), so a
# `grid2d` launch would silently leave every layer past `k = 1` stale rather than error.
# `_driving_rate_bp!` is the column-grid port: `driving_x` has no `z` dependence (broadcast
# read at `k = 1`, the same convention `dotvel!` itself uses), so every layer of the 3D
# scratch gets the identical `driving_x/ρ` value — wasteful relative to a genuinely 2D
# reduction, but correct, and consistent with `solver.residual_x`/`residual_y` already being
# column-shaped for BP regardless of what fills them.

@kernel inbounds = true function _driving_rate_bp!(rate_x, rate_y, driving_x, driving_y, ρ,
                                                    mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    Z = zero(eltype(rate_x))
    rate_x[I...] = node_active(mask, NODE_ACX, i, j) ? driving_x[i, j, 1] / ρ : Z
    rate_y[I...] = node_active(mask, NODE_ACY, i, j) ? driving_y[i, j, 1] / ρ : Z
end

function _convergence_scale(::ScaledResidual, mech::MechanicState, c::Constants,
                            solver::PseudoTransientSolver, rt::Runtime,
                            ::MomentumBalance3D, mask::AbstractIceMask)
    T = eltype(solver.residual_x)
    rt.launch(rt.arch, rt.grid,
              _driving_rate_bp! => (solver.residual_x, solver.residual_y,
                                    mech.stress.driving_x, mech.stress.driving_y,
                                    convert(T, c.density_ice), mask, rt.grid))
    scale = max(maximum(abs, asarray(solver.residual_x)),
                maximum(abs, asarray(solver.residual_y)))
    # A domain with no driving stress at all is already in balance; falling back to 1 makes
    # `err` the raw residual rather than `Inf` — same fallback as the 2D method.
    return scale > zero(scale) ? scale : one(scale)
end

###############################################################
# Vertical-implicit line relaxation (Phase 2)
###############################################################
#
# `roadmaps/blatter-pattyn.md`, Phase 2. What changes is the *step*, not the residual: nothing
# above this point moves, `dotvel!` still writes the fully explicit `r̃(u^k)`, and
# `ImplicitVertical` is read as a **preconditioner** on the velocity update alone.
#
# Explicit (Phase 1, `pseudo_vel!`) the update is `Δu = θ_v Δτ dv`, i.e. `Δu = θ_v·scale·M⁻¹dv`
# with `M = diag(Λ)` the full Gershgorin diagonal and `Δτ = scale/Λ`. Writing `Ã_v` for the
# vertical part of the residual operator (`-∂z(σxz)/ρ̃`, bed BC included, positive
# semi-definite) and `Λ_horiz/ρ̃` for the membrane row sum alone, this phase keeps that shape
# and replaces the preconditioner by
#
#   M = diag(Λ_horiz/ρ̃) + Ã_v,      M Δu = θ_v · scale · dv.
#
# The accumulator `dv` is untouched, which is what makes this a preconditioner and not a
# different discretization — `dotvel!` still writes the true, fully explicit `r̃(u^k)`, and the
# residual history the damping carries stays the true one.
#
# Equivalently and more familiarly: substituting `ℒ_v(u^{k+1})` for `ℒ_v(u^k)` in the update
# gives `(I - θ_v Δτ ℒ_v)Δu = θ_v Δτ dv`, since `ℒ_v` is linear in `u` (the bed flux is exactly
# `β u[1]`). The two coincide at `θ_v·scale = 1` and differ only in how strongly the vertical
# operator is inverted. **The preconditioner form above is the one implemented**, for two
# reasons: `λ_max(M⁻¹Ã) ≤ 1` — what `AutotunedDynamicRelaxation` needs — holds exactly, since
# `M - Ã = diag(Λ_horiz/ρ̃) - Ã_h` is precisely the horizontal Gershgorin defect; and it does
# not degenerate back to an explicit vertical step as `scale` shrinks, which the other form
# does and which would silently give the phase's benefit away whenever the tuner tightens `Δτ`.
#
# Row `k`, divided through by `Λ_horiz/ρ̃` so the diagonal is `1` plus a dimensionless
# stiffness (`τ̂ = 1/Λ_horiz`, cancelling the `1/ρ̃` both sides carry; `θ_v·scale·τ̂_k·ρ̃` is just
# `θ_v Δτ_k`):
#
#   -τ̂ A_k Δu[k-1] + (1 + τ̂(A_k + C_k + D_k)) Δu[k] - τ̂ C_k Δu[k+1] = θ_v Δτ_k dv[k]
#
# with `A_k` the bed-ward interface coupling `_vshear_*(…, k, k)` (absent at `k = 1`), `C_k`
# the surface-ward one `_vshear_*(…, k, k+1)` (absent at `k = nz`, the stress-free surface),
# and `D_k` the bed drag `_vdrag_*`, diagonal-only and present at `k = 1` alone — the same
# three groups `_pseudo_dt_gershgorin_bp!` sums, read off the same helpers so the bound and the
# operator cannot drift apart. All three are pre-mass-scaling, like the helpers themselves, so
# `ρ̃` never appears in this kernel at all.
#
# The diagonal is `≥ 1` and every row is diagonally dominant by construction, so the Thomas
# sweep needs no pivoting and cannot divide by zero.
#
# One work-item per horizontal face, serial in `k`, `rt.launch2d` + an internal `for k` loop —
# the same shape as `_verticalvelocity!` and `_viscosity_integrals!`, for the same reason. No
# column reads another column's velocity, so writing the forward sweep's `d'` into the
# velocity field itself is safe and saves an array; only `c'` needs storage that outlives the
# sweep, and that is what `ImplicitVertical` carries.

@inline function _line_relax_x!(u, u_old, dv, dtau, cp, μ, H, β, θ, drag, nz,
                                mask, dx, dy, grid, grid2d, i, j)
    T = eltype(u)
    Z = zero(T)
    node_active(mask, NODE_ACX, i, j) || return nothing
    Hx = lerp(H, NODE_ACX, grid2d, i, j, 1)
    Hx > Z || return nothing

    # Forward sweep. `b`/`d` never outlive one layer, so only `c'` reaches an array; `u`
    # temporarily holds `d'`, overwritten by the back substitution below. The sub-diagonal is
    # `a_k = -low_k`, hence the `+ low` where the textbook Thomas has `- a`.
    for k in 1:nz
        Λh = _lambda_horiz_x(μ, mask, grid, dx, dy, i, j, k)
        # `Λh == 0` (every neighbouring `µ` masked away) is the degenerate row the explicit
        # scheme handles by `dtau = 0`, i.e. a frozen face; `τ̂ = 0` freezes it here too, and
        # keeps the `1/Λh` from becoming an `Inf` that would poison the whole column.
        τ̂  = Λh > Z ? inv(Λh) : Z
        low = k > 1  ? τ̂ * _vshear_x(μ, mask, grid, Hx, i, j, k, k) : Z
        up  = k < nz ? τ̂ * _vshear_x(μ, mask, grid, Hx, i, j, k, k + 1) : Z
        bed = (k == 1 && drag) ? τ̂ * _vdrag_x(β, grid, grid2d, Hx, i, j) : Z

        b = one(T) + low + up + bed
        d = θ * dtau[i, j, k] * dv[i, j, k]
        w = k > 1 ? b + low * cp[i, j, k - 1] : b
        cp[i, j, k] = -up / w
        u[i, j, k]  = (k > 1 ? d + low * u[i, j, k - 1] : d) / w
    end

    # Back substitution, carrying `Δu[k+1]` in a register so `u` can take its final value
    # (`u_old + Δu`) in the same pass.
    x = u[i, j, nz]
    u[i, j, nz] = u_old[i, j, nz] + x
    for k in (nz - 1):-1:1
        x = u[i, j, k] - cp[i, j, k] * x
        u[i, j, k] = u_old[i, j, k] + x
    end
    return nothing
end

@inline function _line_relax_y!(v, v_old, dv, dtau, cp, μ, H, β, θ, drag, nz,
                                mask, dx, dy, grid, grid2d, i, j)
    T = eltype(v)
    Z = zero(T)
    node_active(mask, NODE_ACY, i, j) || return nothing
    Hy = lerp(H, NODE_ACY, grid2d, i, j, 1)
    Hy > Z || return nothing

    for k in 1:nz
        Λh = _lambda_horiz_y(μ, mask, grid, dx, dy, i, j, k)
        τ̂  = Λh > Z ? inv(Λh) : Z
        low = k > 1  ? τ̂ * _vshear_y(μ, mask, grid, Hy, i, j, k, k) : Z
        up  = k < nz ? τ̂ * _vshear_y(μ, mask, grid, Hy, i, j, k, k + 1) : Z
        bed = (k == 1 && drag) ? τ̂ * _vdrag_y(β, grid, grid2d, Hy, i, j) : Z

        b = one(T) + low + up + bed
        d = θ * dtau[i, j, k] * dv[i, j, k]
        w = k > 1 ? b + low * cp[i, j, k - 1] : b
        cp[i, j, k] = -up / w
        v[i, j, k]  = (k > 1 ? d + low * v[i, j, k - 1] : d) / w
    end

    x = v[i, j, nz]
    v[i, j, nz] = v_old[i, j, nz] + x
    for k in (nz - 1):-1:1
        x = v[i, j, k] - cp[i, j, k] * x
        v[i, j, k] = v_old[i, j, k] + x
    end
    return nothing
end

@kernel inbounds = true function _vertical_line_relax_bp!(ux, uy, ux_old, uy_old, dvx, dvy,
                                                           dtau_x, dtau_y, cpx, cpy, μ, H, β,
                                                           θ, drag, nz, mask, dx, dy,
                                                           grid, grid2d, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    _line_relax_x!(ux, ux_old, dvx, dtau_x, cpx, μ, H, β, θ, drag, nz,
                   mask, dx, dy, grid, grid2d, i, j)
    _line_relax_y!(uy, uy_old, dvy, dtau_y, cpy, μ, H, β, θ, drag, nz,
                   mask, dx, dy, grid, grid2d, i, j)
end

"""
$(TYPEDSIGNATURES)

Advance `ux`/`uy` by one pseudo-transient step, dispatching on
`solver.vertical_treatment`. [`ExplicitVertical`](@ref) is the pair of [`pseudo_vel!`](@ref)
calls the [`MomentumBalance2D`](@ref) loop also makes, unchanged;
[`ImplicitVertical`](@ref) replaces them with a per-column tridiagonal solve of the
vertical-shear operator (see the source note above).

Called from [`pseudo_transient!`](@ref) after the `u → u_old` copy, so `ux_old`/`uy_old` hold
the iterate the increment is measured from.
"""
function _velocity_update!(::ExplicitVertical, solver::PseudoTransientSolver,
                           mech::MechanicState, c::Constants, rt::Runtime,
                           ::AbstractIceMask, ux, uy, ux_old, uy_old, theta_v)
    pseudo_vel!(asarray(ux), asarray(ux_old), asarray(solver.velocity_x_dt),
                asarray(solver.dtau_x), theta_v)
    pseudo_vel!(asarray(uy), asarray(uy_old), asarray(solver.velocity_y_dt),
                asarray(solver.dtau_y), theta_v)
    return nothing
end

function _velocity_update!(vt::ImplicitVertical, solver::PseudoTransientSolver,
                           mech::MechanicState, c::Constants, rt::Runtime,
                           mask::AbstractIceMask, ux, uy, ux_old, uy_old, theta_v)
    T = eltype(solver.dtau_x)
    dx = Δx(rt.grid2d, Center(), 1, 1, 1)
    dy = Δy(rt.grid2d, Center(), 1, 1, 1)
    nz = size(rt.grid, Center())[3]
    rt.launch2d(rt.arch, rt.grid2d,
                _vertical_line_relax_bp! =>
                    (ux, uy, ux_old, uy_old, solver.velocity_x_dt, solver.velocity_y_dt,
                     solver.dtau_x, solver.dtau_y, vt.thomas_x, vt.thomas_y,
                     mech.material.viscosity, mech.topography.thickness,
                     mech.friction.beta_eff, convert(T, theta_v),
                     _drag_in_spectrum(solver.friction_update), nz, mask,
                     convert(T, dx), convert(T, dy), rt.grid, rt.grid2d))
    return nothing
end

# The `M`-inner product `Δuᵀ M Δu`, up to the factor `scale`, for a preconditioner that is no
# longer diagonal. `_sum_du2_over_dtau` evaluates it as `Σ Δu²/Δτ`, which is `Δuᵀ diag(Λ) Δu`
# only because the explicit step *is* `Δu = Δτ dv`. The implicit step is `M Δu = scale·dv` by
# construction (see the source note above), so the same quantity is `Σ Δu·dv` — no vertical
# operator applied a second time, no extra field, and identical to the explicit expression
# wherever both are defined (`θ_v = 1` under `AutotunedDynamicRelaxation`, so
# `Δu·dv = Δu²/Δτ` there). It also needs no `dtau > 0` guard: off-mask both factors are
# exactly zero rather than forming a `0/0`.
@inline _du_dot_dv(u, u_old, dv) = (u - u_old) * dv

_sum_du_dot_dv(solver, ux, uy) =
    mapreduce(_du_dot_dv, +, asarray(ux), asarray(solver.velocity_x_old),
              asarray(solver.velocity_x_dt)) +
    mapreduce(_du_dot_dv, +, asarray(uy), asarray(solver.velocity_y_old),
              asarray(solver.velocity_y_dt))

"""
$(TYPEDSIGNATURES)

[`AbstractVerticalTreatment`](@ref)-aware [`_arm_tuning`](@ref), used by the
[`MomentumBalance3D`](@ref) loop. [`ExplicitVertical`](@ref) delegates to the method the
depth-averaged loop uses, unchanged; [`ImplicitVertical`](@ref) differs in one term, the
Rayleigh quotient's `M`-inner product, which is `Σ Δu·dv` rather than `Σ Δu²/Δτ` once `M`
carries the vertical operator (see `_sum_du_dot_dv`).
"""
_arm_tuning(tu::AbstractPTTuning, ::AbstractVerticalTreatment, state, solver, ux, uy,
            iter::Int) = _arm_tuning(tu, state, solver, ux, uy, iter)

function _arm_tuning(tu::AutotunedDynamicRelaxation, ::ImplicitVertical, state,
                     solver::PseudoTransientSolver, ux, uy, iter::Int)
    iter % tu.cadence == 0 || return state
    return merge(state, (; rayleigh_ur = _sum_du_dot_r(solver, ux, uy),
                           rayleigh_uu = _sum_du_dot_dv(solver, ux, uy),
                           armed = true))
end

"""
$(TYPEDSIGNATURES)

Iterate `mech.velocity.x`/`y` in pseudo-time until the Blatter-Pattyn momentum-balance
residual vanishes — the [`MomentumBalance3D`](@ref) counterpart of
[`pseudo_transient!(..., ::MomentumBalance2D, ...)`](@ref). Same return shape, same three
solver-held strategies, same halo-filling contract; the differences are exactly what
`roadmaps/blatter-pattyn.md` §1 lists:

 - the unknown is `mech.velocity.x`/`y` (`ACX3`/`ACY3`), not the depth-averaged
   `velocity.depthaverage_x`/`y` — BP resolves the full column, so there is nothing to
   reconstruct afterward;
 - the halo refresh ([`bc!`](@ref)) runs on `rt.grid`, the column grid;
 - no DIVA-style depth-integrated-viscosity chain — BP has none to refresh
   (`roadmaps/blatter-pattyn.md`, §2.2), so the loop body is one line shorter;
 - the velocity update goes through [`_velocity_update!`](@ref) rather than
   [`pseudo_vel!`](@ref) directly, so `solver.vertical_treatment` can replace the explicit
   step with a per-column implicit line solve ([`ImplicitVertical`](@ref), Phase 2). Under
   the default [`ExplicitVertical`](@ref) it *is* the same pair of `pseudo_vel!` calls.

Requires a column [`StaggeredGrid`](@ref) (`nz > 1`), checked by
[`_check_momentum_grid`](@ref) — on `nz == 1` BP would silently degenerate to SSA with the
wrong (3D, unaveraged) viscosity rather than raise a clear error.

`strainrate_cap`/`dtau_cap` (both default `Inf`, a no-op) are the numerical safety nets
documented at [`clamp_velocity_gradients!`](@ref) and above
[`_pseudo_dt_gershgorin_bp!`](@ref) respectively — added after real 8 km AIS geometry showed
a masked, ice-free-adjacent column can lose its Gershgorin bound's dominant vertical term
entirely (`roadmaps/blatter-pattyn.md`, Phase 4 "Thin and ice-free columns"). Two
complementary nets, not one: `dtau_cap` bounds the *first* explicit step at such a column,
`strainrate_cap` bounds the membrane-stress feedback a first step that is still too large
would otherwise feed into every following iteration.
"""
function pseudo_transient!(mech::MechanicState, c::Constants, solver::PseudoTransientSolver,
                           rt::Runtime,
                           momentum::MomentumBalance3D,
                           mask::AbstractIceMask = NoMask();
                           strainrate_cap = Inf, dtau_cap = Inf)
    _check_momentum_grid(momentum, rt)

    (; velocity) = mech
    (; abstol, maxiter, printout_every, ncheck, tuning) = solver
    vertical = solver.vertical_treatment

    ux, uy = velocity.x, velocity.y
    ux_old, uy_old = solver.velocity_x_old, solver.velocity_y_old
    resid_x, resid_y = solver.residual_x, solver.residual_y

    # Fields held fixed over the PT iteration.
    drivingstress!(mech, c, rt, momentum, mask)
    state = _tuning_init!(tuning, solver, mech, c, rt, momentum, mask; dtau_cap)
    # After `drivingstress!` (it reads the driving stress) and before the loop (it borrows
    # `residual_x`/`residual_y` as scratch, which `dotvel!` overwrites on iteration 1).
    scale = _convergence_scale(solver.convergence, mech, c, solver, rt, momentum, mask)

    T = eltype(asarray(ux))
    err  = typemax(T)
    iter = 0
    while err > abstol && iter < maxiter
        iter += 1

        pseudo_rate!(mech, c, rt, momentum, solver, mask; gamma = state.gamma, strainrate_cap)
        state = _tune!(tuning, state, solver, mech, c, rt, momentum, mask, ux, uy; dtau_cap)

        # After `_tune!`, which needs the pre-copy `u_old` (see its docstring).
        copyto!(asarray(ux_old), asarray(ux))
        copyto!(asarray(uy_old), asarray(uy))

        _velocity_update!(vertical, solver, mech, c, rt, mask, ux, uy, ux_old, uy_old,
                          state.theta_v)

        bc!(rt.arch, rt.grid, ux => Neumann())
        bc!(rt.arch, rt.grid, uy => Neumann())

        state = _arm_tuning(tuning, vertical, state, solver, ux, uy, iter)

        if iter % ncheck == 0 || iter == maxiter
            err = _pt_error(solver.convergence, solver, ux, uy, scale)
        end
        if iter % printout_every == 0
            println("PT iteration $iter: err = $err, gamma = $(state.gamma)")
        end
    end

    # Diagnostic only (does not gate the loop above): the raw, undamped rate norm at the
    # final iterate, independent of `gamma` — see `dotvel!`'s `resid_x`/`resid_y` note.
    residual = max(maximum(abs, asarray(resid_x)), maximum(abs, asarray(resid_y)))

    return (; iterations = iter, error = err, converged = err <= abstol, residual,
              damping = state.gamma, lambda_min = state.lambda_min)
end

###############################################################
# Mass continuity on the staggered grid
###############################################################
#
# Depth-integrated mass conservation for an ice sheet:
#
#   ∂H/∂t = -∇·q + ṁ,        q = H ū
#
# with H the ice thickness, ū the depth-averaged horizontal velocity and ṁ the net mass
# balance (surface + basal + calving). On the Arakawa C-grid this is a *flux-form*
# discretization: the two flux components live on the cell faces the divergence reads
# (`q_x` at `acx`, `q_y` at `acy`), so
#
#   (∇·q)[i,j] = (q_x[i+1,j] - q_x[i,j])/Δx + (q_y[i,j+1] - q_y[i,j])/Δy
#
# telescopes exactly over any block of cells: whatever leaves one cell enters its
# neighbour, and the total thickness change over a region equals the net flux across its
# boundary plus the mass balance inside it — to machine precision, not to truncation
# order. That property is the reason for the staggering, and it is why the flux fields
# are stored at faces ([`FluxState`](@ref)) rather than re-staggered inside the
# divergence.
#
# What the schemes differ in is only the *reconstruction of H at the face*, since H is a
# cell-centred quantity and the flux is not.

"""
$(TYPEDSIGNATURES)

How the ice thickness is reconstructed on the cell faces when forming the mass flux
`q = H ū`. All implemented schemes are conservative flux-form discretizations of
`∂H/∂t = -∇·q + ṁ`; they differ only in the face value of `H`.

Subtypes: [`NoAdvection`](@ref), [`CenteredAdvection`](@ref), [`UpwindAdvection`](@ref).
"""
abstract type AbstractAdvection end

"""
$(TYPEDSIGNATURES)

No mass transport: the fluxes are zeroed, so `∂H/∂t = ṁ` (a purely local mass balance).
Useful for isolating the mass-balance forcing in tests and for prescribed-geometry runs.
"""
struct NoAdvection <: AbstractAdvection end

"""
$(TYPEDSIGNATURES)

Centred (arithmetic-mean) thickness at the face: `H_face = (H[i-1] + H[i]) / 2`, via
Chmy's `lerp`. Second-order accurate and non-dissipative, but not monotone — steep
margins can produce over/undershoots, and it does not respect the direction of flow.
"""
struct CenteredAdvection <: AbstractAdvection end

"""
$(TYPEDSIGNATURES)

First-order upwind thickness at the face: `H_face` is taken from the cell the ice is
flowing *out of* (`H[i-1]` where the face velocity is positive, `H[i]` where it is
negative). Monotone and positivity-friendly at the cost of numerical diffusion — the
usual choice at ice margins and grounding lines, where [`CenteredAdvection`](@ref)
oscillates.
"""
struct UpwindAdvection <: AbstractAdvection end

"""
$(TYPEDSIGNATURES)

Placeholder: level-set (implicit front-tracking) advection of the ice margin. Not
implemented — no [`mass_flux!`](@ref) method dispatches on it yet.
"""
struct LevelSetAdvection <: AbstractAdvection end

###############################################################
# Face reconstruction of H
###############################################################

# `to` is the destination node class, `dim` the axis being staggered along. Both are
# compile-time constants at the call site, so the branch below is resolved by dispatch,
# not at runtime, and the whole reconstruction inlines into the flux kernel.

@inline _face_thickness(
    ::CenteredAdvection,
    H,
    u_face,
    to,
    dim,
    grid,
    I::Vararg{Integer,3},
) = lerp(H, to, grid, I...)

# `left`/`right` on a Center field with a Vertex destination are Chmy's own names for
# "the cell below / above this face": `left(H, Dim(1), i, j, k) == H[i-1, j, k]` and
# `right(...) == H[i, j, k]`, because Chmy places vertex `i` between centres `i-1` and
# `i`. Using them (rather than literal index arithmetic) keeps the staggering convention
# in one place — Chmy's — instead of duplicating it here.
@inline _face_thickness(
    ::UpwindAdvection,
    H,
    u_face,
    to,
    dim,
    grid,
    I::Vararg{Integer,3},
) = u_face > zero(u_face) ? left(H, dim, I...) : right(H, dim, I...)

###############################################################
# Kernels
###############################################################

@kernel inbounds = true function _mass_flux!(q_x, q_y, H, u, v, scheme, mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    Z = zero(eltype(q_x))

    if node_active(mask, NODE_ACX, i, j)
        u_face = u[I...]
        q_x[I...] =
            _face_thickness(scheme, H, u_face, NODE_ACX, Dim(1), grid, I...) * u_face
    else
        q_x[I...] = Z
    end

    if node_active(mask, NODE_ACY, i, j)
        v_face = v[I...]
        q_y[I...] =
            _face_thickness(scheme, H, v_face, NODE_ACY, Dim(2), grid, I...) * v_face
    else
        q_y[I...] = Z
    end
end

@kernel inbounds = true function _zero_flux!(q_x, q_y, O)
    I = @index(Global, NTuple)
    I = I + O
    q_x[I...] = zero(eltype(q_x))
    q_y[I...] = zero(eltype(q_y))
end

@kernel inbounds = true function _thickness_rate!(dHdt, q_x, q_y, mb, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    dHdt[I...] = _mass_balance(mb, I...) - (∂x(q_x, grid, I...) + ∂y(q_y, grid, I...))
end

# The mass balance is a field in a real run and a scalar in idealized ones; both resolve
# at compile time, so neither costs a branch in the kernel.
@inline _mass_balance(mb::Chmy.AbstractField, I::Vararg{Integer,3}) = mb[I...]
@inline _mass_balance(mb::Number, ::Vararg{Integer,3}) = mb

###############################################################
# Driver functions
###############################################################

"""
$(TYPEDSIGNATURES)

Fill the face-staggered depth-integrated mass flux `q = H ū` from the cell-centred ice
thickness `H` (at `aa`) and the depth-averaged velocity components `u`/`v` (at
`acx`/`acy`), reconstructing `H` at the faces according to `scheme`.

`q_x` must live at `acx` and `q_y` at `acy`, on the *depth-integrated* grid — this is a
2D balance, so it is launched with `rt.launch2d`/`rt.grid2d` (see [`Runtime`](@ref)).

`mask` ([`AbstractIceMask`](@ref), default [`NoMask`](@ref)) is applied **to the fluxes,
which is the only conservative place to apply it.** A face the mask switches off carries no
flux, so no mass crosses it and none is created or destroyed — the divergence still
telescopes exactly. Masking the *tendency* instead would let mass leave a donor cell and
never arrive, silently breaking conservation.

!!! note "Do not mask this by ice presence alone if the margin must advance"
    A face with ice on one side only is exactly the face that carries ice into the empty
    cell, so [`IceMask`](@ref)'s "active if **any** adjoining cell is active" rule keeps
    it open. What a mask here *can* legitimately block is a hard domain constraint
    (`allowed`), which uses the stricter "**all** adjoining cells permitted" rule — ice
    then piles up against the wall rather than being destroyed.
"""
function mass_flux!(
    q_x,
    q_y,
    H,
    u,
    v,
    scheme::AbstractAdvection,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _mass_flux! => (q_x, q_y, H, u, v, scheme, mask, rt.grid2d),
    )
    return nothing
end

function mass_flux!(
    q_x,
    q_y,
    H,
    u,
    v,
    ::NoAdvection,
    rt::Runtime,
    ::AbstractIceMask = NoMask(),
)
    rt.launch2d(rt.arch, rt.grid2d, _zero_flux! => (q_x, q_y))
    return nothing
end

"""
$(TYPEDSIGNATURES)

Compute the ice-thickness tendency `dHdt = -∇·q + ṁ` at `aa` from the face fluxes
`q_x`/`q_y` and the net mass balance `mb` (a `Field` at `aa`, or a scalar).

This is a *rate*, not a step: integrating it in time is the time stepper's job (see
[`Integrator`](@ref)), which is also where any positivity enforcement on `H` belongs —
nothing here prevents a large `Δt` from driving `H` negative.

`bc` is forwarded to Chmy's `Launcher`. Without it only `interior(dHdt)` is meaningful:
the launcher sweeps one halo ring, and at that ring the divergence reads a flux cell the
flux kernel's own sweep did not reach. That is harmless for the usual sequence (step `H`
over the interior, then apply boundary conditions to `H`), but a consumer that reads the
tendency's halo must fill it.

!!! note "Boundary fluxes come from the halo, not from here"
    Whether the domain edge is closed (no outflow) or open follows entirely from the
    halo values of `H` and `ū` when [`mass_flux!`](@ref) ran — the divergence itself makes
    no boundary decision. Mapping ice-sheet boundary conditions onto halo fills is Phase 4
    of `pagos-roadmaps/chmy.md`.
"""
function thickness_rate!(dHdt, q_x, q_y, mb, rt::Runtime; bc = nothing)
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _thickness_rate! => (dHdt, q_x, q_y, mb, rt.grid2d);
        bc,
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

Mass continuity in one call: reconstruct the face fluxes with `scheme`
([`mass_flux!`](@ref)), then take their divergence against the mass balance
([`thickness_rate!`](@ref)). Writes `q_x`, `q_y` and `dHdt`; reads `H`, `u`, `v`, `mb`.

`mask` reaches [`mass_flux!`](@ref) only; the divergence is deliberately **never** masked,
both to keep the scheme conservative and because a cell has to be allowed a tendency before
it holds any ice, or the margin could never advance.
"""
function advect!(
    dHdt,
    q_x,
    q_y,
    H,
    u,
    v,
    mb,
    scheme::AbstractAdvection,
    rt::Runtime,
    mask::AbstractIceMask = NoMask();
    bc = nothing,
)
    mass_flux!(q_x, q_y, H, u, v, scheme, rt, mask)
    thickness_rate!(dHdt, q_x, q_y, mb, rt; bc)
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level [`advect!`](@ref): takes the geometry and mass balance from `topo`, the
depth-averaged velocity from `mech`, stores the fluxes in `mech.flux` and the tendency in
`topo.thickness.ice_dt`.
"""
function advect!(
    topo::TopographicState,
    mech::MechanicState,
    scheme::AbstractAdvection,
    rt::Runtime,
    mask::AbstractIceMask = NoMask();
    bc = nothing,
)
    return advect!(
        topo.thickness.ice_dt,
        mech.flux.x,
        mech.flux.y,
        topo.thickness.ice,
        mech.velocity.depthaverage_x,
        mech.velocity.depthaverage_y,
        topo.massbalance.net,
        scheme,
        rt,
        mask;
        bc,
    )
end

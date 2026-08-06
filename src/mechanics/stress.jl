"""
$(TYPEDSIGNATURES)

Compute the deviatoric stress tensor of `mech` from its strain-rate tensor and the
ice viscosity `mat.eta_ice`, following the constitutive relation `` \\tau_{ij} = 2 \\, \\eta \\, \\dot\\varepsilon_{ij} ``,
and the effective (second-invariant) stress `mech.stress_effective`.

The computation is **independent of the dynamics**: the strain-rate, stress and
viscosity fields are all 3D `(nx, ny, nz)` arrays (depth-averaged solvers simply use
`nz == 1`), so every component is obtained the same way. In particular `` \\tau_{zz} ``
is reconstructed from the traceless property of the deviatoric stress tensor,
`` \\tau_{zz} = -(\\tau_{xx} + \\tau_{yy}) ``, which under incompressibility equals
`` 2\\eta\\dot\\varepsilon_{zz} `` exactly, so reconstructing loses no information. The
vertical-shear components `` \\tau_{xz}, \\tau_{yz} `` follow directly from
`` \\dot\\varepsilon_{xz}, \\dot\\varepsilon_{yz} ``, which are zero in the depth-averaged
case (no resolved vertical shear) and populated in the full-column case.

A **single fused KernelAbstractions kernel** does the whole tensor in one pass on
whichever backend (CPU, GPU) owns the `mech` arrays: it reads the viscosity and the
five independent strain-rate components once, then writes all six stress components
and the effective stress. Because every field is a 3D array of the same shape, the
kernel uses a flat `@index(Global, Linear)` index, so no 2D/3D special-casing is
needed; the single read/write sweep over memory is the bandwidth-optimal layout for
this otherwise memory-bound computation.

This assumes the `strain_rate_d*` fields of `mech` are already populated.
"""
function deviatoric_stress!(mech::MechanicState, mat::MaterialState)
    deviatoric_stress!(
        mech.stress.xx,
        mech.stress.yy,
        mech.stress.zz,
        mech.stress.xy,
        mech.stress.xz,
        mech.stress.yz,
        mech.stress.effective,
        mat.eta_ice,
        mech.strainrate.xx,
        mech.strainrate.yy,
        mech.strainrate.xy,
        mech.strainrate.xz,
        mech.strainrate.yz,
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

Guard: the *collocated* [`deviatoric_stress!`](@ref) applied to a `Field`-based
[`MechanicState`](@ref) would be **silently wrong rather than an error**, so this
combination is rejected outright and points at the staggered method
[`deviatoric_stress!(mech, mat, rt::Runtime)`](@ref) instead.

The fused kernel it would run sweeps a single flat `@index(Global, Linear)` shared by every
tensor component, which is valid exactly while they all have the same shape — true for
plain arrays, false on the C-grid. There the six components sit at four different node
classes with three different lengths (for `nz = 6` on an 8×8 grid: `τ_xx` at `aa` is
8×8×6, `τ_xy` at `ab` is 9×9×6, `τ_xz` at `acx`/z-`Vertex` is 9×8×7), so one linear index
addresses a *different physical point* in each of them. Because the shortest field drives
`ndrange`, the reads stay in bounds and nothing throws — hence the guard rather than
trusting the failure to be visible.
"""
deviatoric_stress!(mech::MechanicState{<:Chmy.AbstractField}, mat::MaterialState) =
    error(
        "this collocated `deviatoric_stress!` cannot be applied to a `Field`-based " *
        "state: its fused kernel assumes every tensor component has the same shape, " *
        "which is false on the C-grid (they sit at four different node classes), and " *
        "it would silently mix locations rather than error. Use the staggered method " *
        "`deviatoric_stress!(mech, mat, rt::Runtime)` instead.",
    )

###############################################################
# Chmy-native, C-grid staggered deviatoric stress
###############################################################
# Strain rates are already at the right node classes (see `mechanics/strainrate.jl`); only
# the viscosity, cell-centred at `aa`, needs interpolating onto the off-diagonal nodes,
# which it does via `hlerp` (see the docstring below for why harmonic). Because `hlerp`
# averages reciprocals, an unmasked ice-free neighbour gives `NaN` — masking the strain rate
# instead doesn't help (`NaN * 0 == NaN`) — so each off-diagonal term is guarded by the
# strict `node_fully_active` rule (every contributing cell icy), not the permissive rule the
# fluxes and gradients use.
@kernel inbounds = true function _deviatoric_stress_staggered!(
    stress,
    sr,
    η,
    mask,
    grid,
    O,
)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    Z = zero(eltype(stress.xx))

    # `aa` needs no interpolation of η, so the permissive rule is enough there: the cell
    # either has ice or it does not.
    if node_active(mask, NODE_AA, i, j)
        twoη = 2 * η[I...]
        txx = twoη * sr.xx[I...]
        tyy = twoη * sr.yy[I...]
        stress.xx[I...] = txx
        stress.yy[I...] = tyy
        stress.zz[I...] = -(txx + tyy)               # traceless deviatoric identity
    else
        stress.xx[I...] = Z
        stress.yy[I...] = Z
        stress.zz[I...] = Z
    end

    txy =
        node_fully_active(mask, NODE_AB, i, j) ?
        2 * hlerp(η, NODE_AB, grid, I...) * sr.xy[I...] : Z
    txz =
        node_fully_active(mask, NODE_ACX_AC, i, j) ?
        2 * hlerp(η, NODE_ACX_AC, grid, I...) * sr.xz[I...] : Z
    tyz =
        node_fully_active(mask, NODE_ACY_AC, i, j) ?
        2 * hlerp(η, NODE_ACY_AC, grid, I...) * sr.yz[I...] : Z
    stress.xy[I...] = txy
    stress.yx[I...] = txy
    stress.xz[I...] = txz
    stress.zx[I...] = txz
    stress.yz[I...] = tyz
    stress.zy[I...] = tyz
end

@kernel inbounds = true function _deviatoric_stress_effective_staggered!(
    stress,
    mask,
    grid,
    O,
)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    if node_active(mask, NODE_AA, i, j)
        txy = lerp(stress.xy, NODE_AA, grid, I...)
        txz = lerp(stress.xz, NODE_AA, grid, I...)
        tyz = lerp(stress.yz, NODE_AA, grid, I...)
        stress.effective[I...] = sqrt(
            (stress.xx[I...]^2 + stress.yy[I...]^2 + stress.zz[I...]^2) / 2 +
            txy^2 +
            txz^2 +
            tyz^2,
        )
    else
        stress.effective[I...] = zero(eltype(stress.effective))
    end
end

"""
$(TYPEDSIGNATURES)

Chmy-native, C-grid staggered [`deviatoric_stress!`](@ref): `τ_ij = 2 η ε̇_ij` with each
component written at its own node class, and the cell-centred viscosity `mat.eta_ice`
interpolated onto the off-diagonal nodes by **harmonic** mean (`hlerp`) — the right average
for a viscosity across which stress, not strain rate, is continuous.

Distinguished from the collocated method by taking a [`Runtime`](@ref). Requires
`mech.strainrate` to be populated (see [`raw_strainrate!`](@ref)). `τ_zz` comes from the
traceless identity `-(τ_xx + τ_yy)`, as in the collocated version.

Runs two launches: the components, then the effective stress at `aa`, which reads
off-diagonal neighbours and so cannot share the first pass. Only
`interior(stress.effective)` is meaningful, for the same reason as
[`raw_strainrate_effective!`](@ref).

!!! warning "A zero viscosity gives `NaN`, not zero"
    `hlerp` averages reciprocals, so `η = 0` in either cell of a pair yields `NaN` rather
    than `0` at the node between them. A freshly allocated [`MaterialState`](@ref) is all
    zeros, so the viscosity must be filled (by a flow law, or explicitly) before this is
    called — with a strictly positive value, which any physical ice viscosity is.
"""
function deviatoric_stress!(
    mech::MechanicState,
    mat::MaterialState,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    rt.launch(
        rt.arch,
        rt.grid,
        _deviatoric_stress_staggered! =>
            (mech.stress, mech.strainrate, mat.eta_ice, mask, rt.grid),
    )
    rt.launch(
        rt.arch,
        rt.grid,
        _deviatoric_stress_effective_staggered! => (mech.stress, mask, rt.grid),
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

Array-level method of [`deviatoric_stress!`](@ref): launches the fused kernel directly
on the stress outputs (`sxx … seff`), the viscosity `η` and the five independent
strain-rate components (`exx, eyy, exy, exz, eyz`). All arrays must share the same shape
and reside on the same backend. `deviatoric_stress!(mech, mat)` is a thin wrapper that
unpacks the corresponding fields of `mech`/`mat` and calls this method.
"""
function deviatoric_stress!(
    sxx,
    syy,
    szz,
    sxy,
    sxz,
    syz,
    seff,
    η,
    exx,
    eyy,
    exy,
    exz,
    eyz,
)
    backend = get_backend(sxx)
    kernel! = _deviatoric_stress_kernel!(backend)
    kernel!(
        sxx,
        syy,
        szz,
        sxy,
        sxz,
        syz,
        seff,
        η,
        exx,
        eyy,
        exy,
        exz,
        eyz;
        ndrange = length(sxx),
    )
    return nothing
end
# @dev TODO this typically does not need to be computed where we don't have any ice and could be easily handled via a mask passed to the kernel. Check performance!

@kernel function _deviatoric_stress_kernel!(
    sxx,
    syy,
    szz,
    sxy,
    sxz,
    syz,
    seff,
    η,
    exx,
    eyy,
    exy,
    exz,
    eyz,
)
    I = @index(Global, Linear)
    @inbounds begin
        twoη = 2 * η[I]
        txx = twoη * exx[I]
        tyy = twoη * eyy[I]
        txy = twoη * exy[I]
        txz = twoη * exz[I]
        tyz = twoη * eyz[I]
        tzz = -(txx + tyy)                      # traceless deviatoric identity
        sxx[I] = txx
        syy[I] = tyy
        szz[I] = tzz
        sxy[I] = txy
        sxz[I] = txz
        syz[I] = tyz
        seff[I] = sqrt((txx^2 + tyy^2 + tzz^2) / 2 + txy^2 + txz^2 + tyz^2)
    end
end

# @dev TODO maybe need to stagger the stresses!
"""
    shearstress!(shear_x, shear_y, strainrate_xx, strainrate_xy, strainrate_yy,
        prealloc, dx, dy)

Compute the shear stress components `shear_x` and `shear_y` from the components
of the scaled strain rate tensor.
"""
function shearstress!(
    shear_x,
    shear_y,
    strainrate_xx,
    strainrate_xy,
    strainrate_yy,
    prealloc,
    dx,
    dy,
)

    ∂x!(prealloc, strainrate_xx, dx)
    shear_x .= prealloc
    ∂y!(prealloc, strainrate_xy, dy)
    shear_x .+= prealloc

    ∂x!(prealloc, strainrate_xy, dx)
    shear_y .= prealloc
    ∂y!(prealloc, strainrate_yy, dy)
    shear_y .+= prealloc
    return nothing
end

"""
$(TYPEDSIGNATURES)

Write the basal-stress components `base_x`, `base_y` from the friction law
`` \\boldsymbol{\\tau}_b = \\beta\\,\\mathbf{v}_b ``. Depth-averaged and
dynamics-independent, so no dispatch on the momentum balance is needed.

Staggering of `β` onto the velocity points is assumed to be handled externally, so the
components are a plain elementwise product. The staggered method below does the staggering
itself.
"""
function basalstress!(base_x, base_y, β, v_x, v_y)
    @. base_x = β * v_x
    @. base_y = β * v_y
    return nothing
end

###############################################################
# Chmy-native, C-grid staggered basal stress
###############################################################
# `β`/`β_eff` live at `aa`, but `τ_b = β v_b` needs `β` on the velocity faces (`acx`/`acy`).
# Resolved by an inline `lerp` (the same choice `drivingstress!` makes for `H`) rather than
# storing `β_acx`/`β_acy` copies, so the face value can never go stale relative to
# `friction.beta_eff`. Arithmetic, not harmonic: unlike the viscosity, `β` has no "two cells
# in series" argument for a harmonic mean, and `lerp` carries no `NaN` risk since `β = 0`
# just means no drag from that side of the face.

@kernel inbounds = true function _basalstress_staggered!(
    base_x,
    base_y,
    β,
    v_x,
    v_y,
    mask,
    grid,
    O,
)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    Z = zero(eltype(base_x))
    base_x[I...] =
        node_active(mask, NODE_ACX, i, j) ? lerp(β, NODE_ACX, grid, I...) * v_x[I...] :
        Z
    base_y[I...] =
        node_active(mask, NODE_ACY, i, j) ? lerp(β, NODE_ACY, grid, I...) * v_y[I...] :
        Z
end

"""
$(TYPEDSIGNATURES)

Chmy-native, C-grid staggered [`basalstress!`](@ref): writes `base_x` at `acx`, `base_y`
at `acy` as `β_face · v_face`, with the cell-centred friction coefficient `β` interpolated
onto the velocity faces by (arithmetic) `lerp` — see the note above for why arithmetic,
not harmonic like the viscosity.

Distinguished from the collocated method by taking a [`Runtime`](@ref). Depth-integrated
throughout, so it runs on `rt.grid2d`.
"""
function basalstress!(
    base_x,
    base_y,
    β,
    v_x,
    v_y,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _basalstress_staggered! => (base_x, base_y, β, v_x, v_y, mask, rt.grid2d),
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level Chmy-native [`basalstress!`](@ref): writes `mech.stress.base_x`/`base_y` from
`mech.friction.beta_eff` and `mech.velocity.base_x`/`base_y`.
"""
basalstress!(mech::MechanicState, rt::Runtime, mask::AbstractIceMask = NoMask()) =
    basalstress!(
        mech.stress.base_x,
        mech.stress.base_y,
        mech.friction.beta_eff,
        mech.velocity.base_x,
        mech.velocity.base_y,
        rt,
        mask,
    )

"""
$(TYPEDSIGNATURES)

Write the driving-stress components `driving_x`, `driving_y` from the surface elevation
`surface`, the thickness `thickness`, the ice density `ρ_ice` and gravity `g`:
`` \\tau_{d} = \\rho_{ice}\\,g\\,H\\,\\nabla s ``, with the surface gradient taken by the
central-difference stencils `∂x!`/`∂y!`. **Independent of the dynamics** — every momentum
balance shares this expression — so no dispatch on the momentum balance is needed.
"""
function drivingstress!(driving_x, driving_y, surface, thickness, ρ_ice, g, dx, dy)
    ∂x!(driving_x, surface, dx)
    ∂y!(driving_y, surface, dy)
    @. driving_x *= ρ_ice * g * thickness
    @. driving_y *= ρ_ice * g * thickness
    return nothing
end

###############################################################
# Chmy-native, C-grid staggered driving stress
###############################################################
# `s` and `H` are cell-centred, but the driving stress must live on the velocity faces it
# forces. Chmy's `∂x`/`∂y` map `aa → acx`/`acy` natively, so the gradient itself needs no
# interpolation (unlike the collocated method's cell-centred central difference); only `H`
# has to be moved, by one `lerp` onto the same face.

@kernel inbounds = true function _surface_gradient!(dsdx, dsdy, s, mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    Z = zero(eltype(dsdx))
    dsdx[I...] = node_active(mask, NODE_ACX, i, j) ? ∂x(s, grid, I...) : Z
    dsdy[I...] = node_active(mask, NODE_ACY, i, j) ? ∂y(s, grid, I...) : Z
end

@kernel inbounds = true function _drivingstress!(τx, τy, s, H, ρg, mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    Z = zero(eltype(τx))
    τx[I...] =
        node_active(mask, NODE_ACX, i, j) ?
        ρg * lerp(H, NODE_ACX, grid, I...) * ∂x(s, grid, I...) : Z
    τy[I...] =
        node_active(mask, NODE_ACY, i, j) ?
        ρg * lerp(H, NODE_ACY, grid, I...) * ∂y(s, grid, I...) : Z
end

"""
$(TYPEDSIGNATURES)

Fill the surface-elevation gradients `dsdx` (at `acx`) and `dsdy` (at `acy`) from the
cell-centred surface elevation `s`, using Chmy's staggered `∂x`/`∂y` — which map
`aa → acx`/`acy` natively, so no interpolation enters.

Independent of [`drivingstress!`](@ref) by design: the driving stress recomputes the
gradient inline rather than reading these fields, so it stays a single memory sweep and
carries no ordering dependency on this function. Call this one when the gradients
themselves are wanted (diagnostics, an SIA diffusivity), not as a prerequisite.
"""
function surface_gradient!(dsdx, dsdy, s, rt::Runtime, mask::AbstractIceMask = NoMask())
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _surface_gradient! => (dsdx, dsdy, s, mask, rt.grid2d),
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level [`surface_gradient!`](@ref): writes `topo.elevation.surface_dx`/`surface_dy`
from `topo.elevation.surface`.
"""
surface_gradient!(
    topo::TopographicState,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
) = surface_gradient!(
    topo.elevation.surface_dx,
    topo.elevation.surface_dy,
    topo.elevation.surface,
    rt,
    mask,
)

"""
$(TYPEDSIGNATURES)

Chmy-native, C-grid staggered [`drivingstress!`](@ref): writes `τx` at `acx` and `τy` at
`acy` from the cell-centred surface elevation `s` and thickness `H`, as
`ρ_ice · g · H_face · ∂s/∂x`, with `H` interpolated onto the face by `lerp` and the
gradient taken by Chmy's `∂x`/`∂y` (`aa → acx`/`acy`, so it needs no interpolation).

Distinguished from the collocated method by taking a [`Runtime`](@ref) instead of `dx`,
`dy`. Everything here is depth-integrated, so it runs on `rt.grid2d`.

!!! warning "Sign convention: this stores `+ρgH∇s`, not the driving stress itself"
    The physical driving stress is `τ_d = -ρ_ice g H ∇s`. What is stored — matching the
    collocated method it replaces, and the pseudo-transient residual that consumes it
    (`dotvel!` computes `(∇·τ - τ_b - driving) / (ρH)`, i.e. it *subtracts* this field) —
    is `+ρ_ice g H ∇s`. Changing the sign here silently reverses the flow direction of
    every momentum balance.

`ρ_ice * g` is converted to the output eltype before entering the kernel, so a Float32
pipeline stays in Float32 (see `roadmaps/chmy.md`, Float32 discipline).
"""
function drivingstress!(
    τx,
    τy,
    s,
    H,
    ρ_ice,
    g,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
)
    ρg = convert(eltype(τx), ρ_ice * g)
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _drivingstress! => (τx, τy, s, H, ρg, mask, rt.grid2d),
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level [`drivingstress!`](@ref): reads the geometry from `mech.topography` (the
mechanics component's own surface/thickness copies) and writes
`mech.stress.driving_x`/`driving_y`.
"""
drivingstress!(
    mech::MechanicState,
    c::Constants,
    rt::Runtime,
    mask::AbstractIceMask = NoMask(),
) = drivingstress!(
    mech.stress.driving_x,
    mech.stress.driving_y,
    mech.topography.surface,
    mech.topography.thickness,
    c.density_ice,
    c.gravity,
    rt,
    mask,
)

###############################################################
# Chmy-native, C-grid staggered Blatter-Pattyn (un-integrated) driving stress
###############################################################
# `ρg ∂s/∂x`, not `ρgH ∂s/∂x` — see `roadmaps/blatter-pattyn.md` §1 for why no thickness
# enters. Written once on `grid2d` and read at `k = 1` inside the 3D [`dotvel!`](@ref) sweep,
# the same broadcast convention `_velocity_gradients!` uses for `H`. No `H` to stagger here,
# so no `lerp` appears at all.

@kernel inbounds = true function _drivingstress_bp!(τx, τy, s, ρg, mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    Z = zero(eltype(τx))
    τx[I...] = node_active(mask, NODE_ACX, i, j) ? ρg * ∂x(s, grid, I...) : Z
    τy[I...] = node_active(mask, NODE_ACY, i, j) ? ρg * ∂y(s, grid, I...) : Z
end

"""
$(TYPEDSIGNATURES)

Chmy-native, C-grid staggered Blatter-Pattyn driving stress: writes `τx` at `acx` and `τy`
at `acy` from the cell-centred surface elevation `s`, as `ρ_ice · g · ∂s/∂x` — no thickness,
unlike [`drivingstress!(::MomentumBalance2D)`](@ref), since the Blatter-Pattyn residual is
per unit volume. Runs on `rt.grid2d`: the surface slope has no `z` dependence, so this is
written once and read broadcast down the column (see the source note above), not copied into
a 3D field.

Same sign convention as the depth-integrated method: stores `+ρgH∇s` scaled *without* `H`,
i.e. `+ρg∇s`, matching what [`dotvel!`](@ref) subtracts.
"""
function drivingstress!(
    τx,
    τy,
    s,
    ρ_ice,
    g,
    rt::Runtime,
    ::MomentumBalance3D,
    mask::AbstractIceMask = NoMask(),
)
    ρg = convert(eltype(τx), ρ_ice * g)
    rt.launch2d(
        rt.arch,
        rt.grid2d,
        _drivingstress_bp! => (τx, τy, s, ρg, mask, rt.grid2d),
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level [`drivingstress!`](@ref) for the Blatter-Pattyn momentum balance: reads the
geometry from `mech.topography` and writes `mech.stress.driving_x`/`driving_y` — the same
fields the depth-integrated method writes, since both are genuinely 2D.
"""
drivingstress!(
    mech::MechanicState,
    c::Constants,
    rt::Runtime,
    momentum::MomentumBalance3D,
    mask::AbstractIceMask = NoMask(),
) = drivingstress!(
    mech.stress.driving_x,
    mech.stress.driving_y,
    mech.topography.surface,
    c.density_ice,
    c.gravity,
    rt,
    momentum,
    mask,
)

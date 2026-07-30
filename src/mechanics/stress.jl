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
        mech.stress.xx, mech.stress.yy, mech.stress.zz,
        mech.stress.xy, mech.stress.xz, mech.stress.yz, mech.stress.effective,
        mat.eta_ice,
        mech.strainrate.xx, mech.strainrate.yy,
        mech.strainrate.xy, mech.strainrate.xz, mech.strainrate.yz,
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
    error("this collocated `deviatoric_stress!` cannot be applied to a `Field`-based " *
          "state: its fused kernel assumes every tensor component has the same shape, " *
          "which is false on the C-grid (they sit at four different node classes), and " *
          "it would silently mix locations rather than error. Use the staggered method " *
          "`deviatoric_stress!(mech, mat, rt::Runtime)` instead.")

###############################################################
# Chmy-native, C-grid staggered deviatoric stress
###############################################################
#
# `τ_ij = 2 η ε̇_ij`, one node class at a time. The strain rates are already at the right
# places (see `mechanics/strainrate.jl`); what has to move is the viscosity, which lives at
# `aa` and must scale off-diagonal components that do not. That interpolation is
# **harmonic** (`hlerp`), not arithmetic: across a viscosity contrast the physically
# conserved quantity is the stress, so the effective viscosity of two cells in series is
# their harmonic mean. An arithmetic mean lets a stiff cell dominate its soft neighbour and
# over-stiffens the margin — the standard choice in staggered viscous solvers, and the
# reason Chmy provides `hlerp` at all.

# The mask is load-bearing here, not a convenience: `η` is zero where there is no ice, and
# `hlerp` averages reciprocals, so an unmasked off-diagonal component is `NaN` on every
# ice-free node. Masking the *strain rate* does not help — `NaN * 0 == NaN` — so the guard
# has to sit around the viscosity term itself, and it has to skip the evaluation rather than
# discard its result.
#
# It also has to be the **strict** rule (`node_fully_active`, every contributing cell icy),
# not the permissive one the fluxes and gradients use. A corner with even one ice-free
# neighbour still gives `hlerp` a zero to invert. And strict is the physically right answer
# as well as the finite one: the harmonic mean of `(η, 0)` tends to `0`, i.e. a cell of no
# viscosity transmits no stress, so a zero there is the limit rather than a fudge. Callers
# pass one mask and this kernel tightens it itself — the alternative, a separate strict mask
# per call site, is a footgun that reintroduces the `NaN` the moment it is forgotten.
@kernel inbounds = true function _deviatoric_stress_staggered!(stress, sr, η, mask, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    Z = zero(eltype(stress.xx))

    # `aa` needs no interpolation of η, so the permissive rule is enough there: the cell
    # either has ice or it does not.
    if node_active(mask, NODE_AA, i, j)
        twoη = 2 * η[I...]
        txx  = twoη * sr.xx[I...]
        tyy  = twoη * sr.yy[I...]
        stress.xx[I...] = txx
        stress.yy[I...] = tyy
        stress.zz[I...] = -(txx + tyy)               # traceless deviatoric identity
    else
        stress.xx[I...] = Z; stress.yy[I...] = Z; stress.zz[I...] = Z
    end

    txy = node_fully_active(mask, NODE_AB, i, j) ?
          2 * hlerp(η, NODE_AB, grid, I...) * sr.xy[I...] : Z
    txz = node_fully_active(mask, NODE_ACX_AC, i, j) ?
          2 * hlerp(η, NODE_ACX_AC, grid, I...) * sr.xz[I...] : Z
    tyz = node_fully_active(mask, NODE_ACY_AC, i, j) ?
          2 * hlerp(η, NODE_ACY_AC, grid, I...) * sr.yz[I...] : Z
    stress.xy[I...] = txy; stress.yx[I...] = txy
    stress.xz[I...] = txz; stress.zx[I...] = txz
    stress.yz[I...] = tyz; stress.zy[I...] = tyz
end

@kernel inbounds = true function _deviatoric_stress_effective_staggered!(stress, mask,
                                                                        grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, _ = I
    if node_active(mask, NODE_AA, i, j)
        txy = lerp(stress.xy, NODE_AA, grid, I...)
        txz = lerp(stress.xz, NODE_AA, grid, I...)
        tyz = lerp(stress.yz, NODE_AA, grid, I...)
        stress.effective[I...] = sqrt((stress.xx[I...]^2 + stress.yy[I...]^2 +
                                       stress.zz[I...]^2) / 2 + txy^2 + txz^2 + tyz^2)
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
function deviatoric_stress!(mech::MechanicState, mat::MaterialState, rt::Runtime,
                            mask::AbstractIceMask = NoMask())
    rt.launch(rt.arch, rt.grid,
              _deviatoric_stress_staggered! =>
                  (mech.stress, mech.strainrate, mat.eta_ice, mask, rt.grid))
    rt.launch(rt.arch, rt.grid,
              _deviatoric_stress_effective_staggered! => (mech.stress, mask, rt.grid))
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
function deviatoric_stress!(sxx, syy, szz, sxy, sxz, syz, seff, η, exx, eyy, exy, exz, eyz)
    backend = get_backend(sxx)
    kernel! = _deviatoric_stress_kernel!(backend)
    kernel!(sxx, syy, szz, sxy, sxz, syz, seff, η, exx, eyy, exy, exz, eyz;
        ndrange = length(sxx))
    return nothing
end
# TODO this typically does not need to be computed where we don't have any ice and could be easily handled via a mask passed to the kernel. Check performance!

@kernel function _deviatoric_stress_kernel!(
    sxx, syy, szz, sxy, sxz, syz, seff,
    η, exx, eyy, exy, exz, eyz,
)
    I = @index(Global, Linear)
    @inbounds begin
        twoη = 2 * η[I]
        txx  = twoη * exx[I]
        tyy  = twoη * eyy[I]
        txy  = twoη * exy[I]
        txz  = twoη * exz[I]
        tyz  = twoη * eyz[I]
        tzz  = -(txx + tyy)                      # traceless deviatoric identity
        sxx[I] = txx; syy[I] = tyy; szz[I] = tzz
        sxy[I] = txy; sxz[I] = txz; syz[I] = tyz
        seff[I] = sqrt((txx^2 + tyy^2 + tzz^2) / 2 + txy^2 + txz^2 + tyz^2)
    end
end

# TODO maybe need to stagger the stresses!
"""
    shearstress!(shear_x, shear_y, strainrate_xx, strainrate_xy, strainrate_yy,
        prealloc, dx, dy)

Compute the shear stress components `shear_x` and `shear_y` from the components
of the scaled strain rate tensor.
"""
function shearstress!(shear_x, shear_y, strainrate_xx, strainrate_xy, strainrate_yy,
    prealloc, dx, dy)

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

Compute the basal-stress components `stress.base_x` and `stress.base_y` of `m` in place.

The basal stress follows the friction law `` \\boldsymbol{\\tau}_b = \\beta\\,\\mathbf{v}_b ``,
i.e. the effective basal friction coefficient `friction.beta_eff` times the basal velocity
`velocity.x_base`, `velocity.y_base`. Like the driving stress it is a depth-averaged 2D
field and dynamics-independent, so no dispatch on the momentum balance is needed.

Staggering of `beta_eff` onto the velocity points is assumed to be handled externally, so
the components are formed by a plain elementwise product.
"""
function basalstress!(m::Mechanics)
    (; state) = m
    (; stress, friction, velocity) = state
    return basalstress!(
        stress.base_x, stress.base_y,
        friction.beta_eff, velocity.x_base, velocity.y_base,
    )
end

"""
$(TYPEDSIGNATURES)

Array-level method of [`basalstress!`](@ref): writes the basal-stress components `base_x`,
`base_y` as the elementwise product of the effective basal friction coefficient `β` and the
basal velocity components `v_x`, `v_y`.
"""
function basalstress!(base_x, base_y, β, v_x, v_y)
    @. base_x = β * v_x
    @. base_y = β * v_y
    return nothing
end
# TODO: this should be staggered!

"""
$(TYPEDSIGNATURES)

Compute the driving-stress components `stress.driving_x` and `stress.driving_y` of `m`
in place.

The driving stress is `` \\tau_{d} = \\rho_{ice}\\,g\\,H\\,\\nabla s ``, where `` s `` is the
ice surface elevation, `` H `` the thickness and `` \\rho_{ice}, g `` are taken from the
physical `Constants`. It is **independent of the dynamics**: every momentum balance shares
this expression, so the depth-averaged 2D `driving_*` fields are obtained the same way and
no dispatch on the momentum balance is needed.
"""
function drivingstress!(m::Mechanics, c::Constants)
    (; state, grid) = m
    (; stress, topography) = state
    (; dx, dy) = grid
    return drivingstress!(
        stress.driving_x, stress.driving_y,
        topography.surface, topography.thickness,
        c.density_ice, c.gravity, dx, dy,
    )
end

"""
$(TYPEDSIGNATURES)

Array-level method of [`drivingstress!`](@ref): writes the driving-stress components
`driving_x`, `driving_y` from the surface elevation `surface`, the thickness `thickness`,
the ice density `ρ_ice`, the gravitational acceleration `g` and the grid spacings `dx`,
`dy`. The surface gradient is taken with the central-difference stencils `∂x!`/`∂y!`.
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
#
# The driving stress is the field the C-grid exists for. `s` and `H` are cell-centred
# (`aa`), the velocity components live on the faces (`acx`/`acy`), and so must the driving
# stress that forces them — the collocated methods above compute `∂s/∂x` with a central
# difference *at the cell centre*, which needs a half-cell interpolation before it can
# force a face velocity. Staggered, no interpolation is needed for the gradient at all:
# Chmy's `∂x` maps `aa → acx` natively, and the difference it takes is exactly the
# two-point difference across that face. Only `H` has to be moved, by one `lerp` onto the
# same face.

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
    τx[I...] = node_active(mask, NODE_ACX, i, j) ?
               ρg * lerp(H, NODE_ACX, grid, I...) * ∂x(s, grid, I...) : Z
    τy[I...] = node_active(mask, NODE_ACY, i, j) ?
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
function surface_gradient!(dsdx, dsdy, s, rt::Runtime,
                           mask::AbstractIceMask = NoMask())
    rt.launch2d(rt.arch, rt.grid2d,
                _surface_gradient! => (dsdx, dsdy, s, mask, rt.grid2d))
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level [`surface_gradient!`](@ref): writes `topo.elevation.surface_dx`/`surface_dy`
from `topo.elevation.surface`.
"""
surface_gradient!(topo::TopographicState, rt::Runtime,
                 mask::AbstractIceMask = NoMask()) =
    surface_gradient!(topo.elevation.surface_dx, topo.elevation.surface_dy,
                      topo.elevation.surface, rt, mask)

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
function drivingstress!(τx, τy, s, H, ρ_ice, g, rt::Runtime,
                        mask::AbstractIceMask = NoMask())
    ρg = convert(eltype(τx), ρ_ice * g)
    rt.launch2d(rt.arch, rt.grid2d,
                _drivingstress! => (τx, τy, s, H, ρg, mask, rt.grid2d))
    return nothing
end

"""
$(TYPEDSIGNATURES)

State-level Chmy-native [`drivingstress!`](@ref): reads the geometry from
`mech.topography` (the mechanics component's own surface/thickness copies, as the
collocated `Mechanics` method does) and writes `mech.stress.driving_x`/`driving_y`.
"""
drivingstress!(mech::MechanicState, c::Constants, rt::Runtime,
               mask::AbstractIceMask = NoMask()) =
    drivingstress!(mech.stress.driving_x, mech.stress.driving_y,
                   mech.topography.surface, mech.topography.thickness,
                   c.density_ice, c.gravity, rt, mask)
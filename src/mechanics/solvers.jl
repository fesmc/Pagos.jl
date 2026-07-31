###############################################################
# Solvers
##############################################################

"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the dynamics solver via [`velocity`](@ref).

# Available subtypes:
"""
abstract type AbstractMomentumSolver end

# TODO
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via energy minimization.
"""
struct OptimMomentumSolver <: AbstractMomentumSolver
end

# TODO
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via an iterative linear solver (e.g., CG, GMRES).
"""
struct IterativeMomentumSolver <: AbstractMomentumSolver
end

"""
$(TYPEDSIGNATURES)

TODO: Solve the ice dynamics via wavelet methods.
"""
struct WaveletMomentumSolver <: AbstractMomentumSolver
end

# TODO
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via convolutional neural network (IGM style).
"""
struct ConvolutionalMomentumSolver <: AbstractMomentumSolver
end

# TODO
"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via a transient solver (e.g., explicit time-stepping).
"""
struct TransientMomentumSolver <: AbstractMomentumSolver
end

"""
$(TYPEDSIGNATURES)

An abstract type to multiple-dispatch whether and how [`PseudoTransientSolver`](@ref)
updates `material.viscosity_depthaveraged` from the current velocity iterate during the
PT loop, following the codebase's "dispatch, not `if`/`else`" convention
(`roadmaps/pagos.md`). [`NoViscosityContinuation`](@ref) (the default) makes the solver
treat viscosity as a fixed input, exactly as every solver behaved before this type existed.
[`GlenViscosityContinuation`](@ref) enables Sandip et al. (2024) Eq. 8's log-space
continuation.

!!! note "Standalone, not integrated with `AbstractFlowLaw`/`AbstractCreep`"
    `GlenViscosityContinuation` implements Glen's law directly from the strain-rate
    invariant (Sandip Eq. 3), independent of `src/material/flow_law.jl`/`creep.jl`. Those
    are built around a *stress*-driven creep formulation (`creep(σ_e, law)`), which would
    need an inner per-cell nonlinear solve to use here (`σ_e = 2ηε̇_e` depends on the very
    `η` being solved for) — not what Sandip's direct, closed-form update does. This is
    flagged, not resolved, in `roadmaps/PT-autotune.md`: the codebase now has two
    parallel notions of "Glen's law" (this one and the stress-driven one), plus two
    *dead* strain-rate-driven copies (`src/legacy/flow_law.jl`,
    `src/thermodynamics/viscosity.jl`, neither `include`d) that were the closer relatives
    of this one and went unused instead of being reused. Making this coherent — one flow-
    law abstraction, or an explicit, documented reason for two — is future work.
"""
abstract type AbstractViscosityContinuation end

"""
$(TYPEDSIGNATURES)

No viscosity continuation: `material.viscosity_depthaveraged` is a fixed input for the
whole PT solve, untouched by the solver. The default for [`PseudoTransientSolver`](@ref),
and the only behaviour prior to [`AbstractViscosityContinuation`](@ref) existing.
"""
struct NoViscosityContinuation <: AbstractViscosityContinuation end

"""
$(TYPEDSIGNATURES)

Glen-law viscosity continuation (Sandip et al. 2024, Eq. 3 and 8): at every PT iteration,
`material.viscosity_depthaveraged` is updated from the current effective strain rate
`ε̇_e` (Sandip Eq. 3, `η = A^(-1/n)·(ε̇_e² + ε̇₀²)^((1-n)/(2n))/2`, regularized by the
strain-rate floor `strainrate_reg = ε̇₀`) and relaxed toward its own previous value in
log-space (Eq. 8, `η_new = exp(θ_μ·log(η_raw) + (1-θ_μ)·log(η_old))`) rather than applied
directly — direct application is what Sandip found could diverge the PT iteration before
the nonlinear viscosity has had time to relax, `θ_μ` controlling how much time.

# Fields:
 - `n_glen`: Glen's flow-law exponent (default `3`, matching [`GlenNyeCreep`](@ref)'s).
 - `theta_mu`: log-space relaxation weight, `0 < theta_mu ≤ 1` (`1` disables relaxation:
   `η_new` is `η_raw` outright, no memory of `η_old`).
 - `strainrate_reg`: strain-rate floor `ε̇₀` preventing the `ε̇_e = 0` singularity
   (`η → ∞` as `ε̇_e → 0` in Glen's law); must be chosen relative to the problem's own
   strain-rate scale, there is no unit-independent default that means the same thing
   everywhere — the constructor requires it rather than guessing.

Reads `material.rate_factor_depthaveraged` as the (prescribed, not thermally coupled) rate
factor `A` — see that field's docstring.
"""
struct GlenViscosityContinuation{T<:AbstractFloat} <: AbstractViscosityContinuation
    n_glen::T
    theta_mu::T
    strainrate_reg::T
end

function GlenViscosityContinuation(T::Type{<:AbstractFloat} = Float64;
    n_glen = 3, theta_mu = 0.1, strainrate_reg,
)
    return GlenViscosityContinuation{T}(T(n_glen), T(theta_mu), T(strainrate_reg))
end

"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via a pseudo-transient (PT) solver following Sandip et al. (2024).
The velocity field is relaxed in pseudo-time until the momentum balance is satisfied,
which only requires local (stencil) operations. All work arrays live on the backend of
the grid the solver is constructed from, so the same code runs on CPU and GPU.

Construct from a grid, overriding parameters selectively via keyword arguments:

```julia
solver = PseudoTransientSolver(grid)
solver = PseudoTransientSolver(grid; maxiter = 500, abstol = 1e-10)
```

# Fields:
 - `ndim1`, `ndim2`, `ndim3`: numerical-dimensionality constants of the PT time step
   (1D, 2D, 3D stencils).
 - `muB`: bulk-to-shear viscosity ratio entering the PT time step.
 - `theta_v`: relaxation weight of the velocity update.
 - `gamma`: damping coefficient of the pseudo-transient rate (Sandip et al. 2024,
   Eq. 12–14; Frankel 1950). The rate is accumulated as
   `dv_new = (1 - gamma) * dv_old + r(u)`, so `gamma = 1` discards all memory and
   recovers the plain (undamped, first-order) PT iteration — the default, and the only
   value exercised by the collocated solver. `0 < gamma < 1` turns the iteration into a
   damped wave, which is what buys sub-quadratic iteration-count scaling; picking a
   good value currently needs the same manual tuning as `theta_v`. Automatically
   selecting it from spectral estimates is the Duretz et al. (2026) autotuning work
   (`roadmaps/PT-autotune.md`), not yet implemented.
 - `abstol`: convergence tolerance on the max-norm velocity change per iteration.
 - `maxiter`: maximum number of PT iterations.
 - `ncheck`: check convergence every `ncheck` iterations. The check is the only
   operation that forces the host to wait for the device (see
   [Asynchronous kernel launches](@ref)), so on GPU a larger value (10–50) keeps
   the launch pipeline busy at the cost of up to `ncheck - 1` extra iterations.
 - `printout_every`: print a convergence monitor every `printout_every` iterations
   (silent by default).
 - `dtau_scaling`: safety scaling of the PT time step.
 - `velocity_x_old`, `velocity_y_old`: previous velocity iterate.
 - `velocity_x_dt`, `velocity_y_dt`: pseudo-transient velocity rate (the damping
   accumulator described under `gamma`).
 - `dtau_x`, `dtau_y`: pseudo-time step, held at the same node class as the velocity
   component it advances. Only the Chmy-native, C-grid staggered
   [`pseudo_transient!`](@ref) fills these per grid point (from the local viscosity, via
   [`pseudo_dt!`](@ref)); the collocated solver keeps its single global scalar (from
   [`pseudo_dt`](@ref)) and never reads these fields, which is why they may look unused
   in that code path.
 - `residual_x`, `residual_y`: the raw (undamped) pseudo-transient rate `r(u)/(ρH)`,
   written unconditionally by [`dotvel!`](@ref) regardless of `gamma` — a genuinely
   separate buffer from `velocity_x_dt`/`velocity_y_dt` because damping's whole point is
   for the accumulator to *not* equal the raw rate, so the raw rate needs its own home to
   remain readable. The Chmy-native [`pseudo_transient!`](@ref) reduces these into the
   `residual` diagnostic it returns; unused by the collocated solver, like `dtau_x`/`dtau_y`.
 - `viscosity_continuation`: an [`AbstractViscosityContinuation`](@ref), default
   [`NoViscosityContinuation`](@ref). Only read by the Chmy-native, C-grid staggered
   [`pseudo_transient!`](@ref); the collocated solver never updates viscosity regardless
   of this field's value.

!!! note "Two type parameters for the work arrays, not one"
    On a [`StaggeredGrid`](@ref) the x- and y-velocity work arrays are Chmy `Field`s at
    different locations (`acx`/`acy`), which are different concrete types (location is
    part of a `Field`'s type). A single `M` shared by all four fields — as a plain-array
    `RegularGrid` state can get away with, since `Matrix{T}` carries no location — would
    reject that combination outright, hence the `MX`/`MY` split below. They collapse to
    the same type on a `RegularGrid`, so that path is unaffected.
"""
struct PseudoTransientSolver{T<:AbstractFloat, MX, MY, VC<:AbstractViscosityContinuation} <: AbstractMomentumSolver
    ndim1::T
    ndim2::T
    ndim3::T
    muB::T
    theta_v::T
    gamma::T
    abstol::T
    maxiter::Int
    ncheck::Int
    printout_every::Int
    dtau_scaling::T
    velocity_x_old::MX
    velocity_y_old::MY
    velocity_x_dt::MX
    velocity_y_dt::MY
    dtau_x::MX
    dtau_y::MY
    residual_x::MX
    residual_y::MY
    viscosity_continuation::VC
end
Adapt.@adapt_structure PseudoTransientSolver

function PseudoTransientSolver(grid::RegularGrid;
    ndim1 = 2.1,
    ndim2 = 4.1,
    ndim3 = 6.1,
    muB = 1e2,
    theta_v = 0.6,
    gamma = 1,
    abstol = 1e-8,
    maxiter = 100,
    ncheck = 1,
    printout_every = typemax(Int),
    dtau_scaling = 1,
    viscosity_continuation::AbstractViscosityContinuation = NoViscosityContinuation(),
)
    T = eltype(grid.x)
    backend = get_backend(grid.x)
    w() = KernelAbstractions.zeros(backend, T, grid.nx, grid.ny)
    return PseudoTransientSolver(
        T(ndim1), T(ndim2), T(ndim3), T(muB),
        T(theta_v), T(gamma), T(abstol), Int(maxiter), Int(ncheck),
        Int(printout_every), T(dtau_scaling), w(), w(), w(), w(), w(), w(), w(), w(),
        viscosity_continuation,
    )
end

"""
$(TYPEDSIGNATURES)

Build a [`PseudoTransientSolver`](@ref) of Chmy `Field`s for the Chmy-native, C-grid
staggered [`pseudo_transient!`](@ref). The velocity-iterate/rate work arrays
(`velocity_x_old`, `velocity_y_old`, `velocity_x_dt`, `velocity_y_dt`) are placed at
`acx`/`acy` on `grid.grid2d`, matching the velocity components they mirror
(`mech.velocity.x`/`y`).

Requires a depth-averaged grid (`grid.grid2d === grid.grid`, i.e. `nz == 1`), checked at
construction rather than left to fail at first solve: DIVA's vertical-shear integral is
Phase 3 future work (`roadmaps/chmy.md`), so today's Chmy-native `pseudo_transient!` only
supports the depth-averaged case, exactly as the collocated method does in practice
(it reads `view(velocity.x, :, :, 1)`).
"""
function PseudoTransientSolver(grid::StaggeredGrid;
    ndim1 = 2.1,
    ndim2 = 4.1,
    ndim3 = 6.1,
    muB = 1e2,
    theta_v = 0.6,
    gamma = 1,
    abstol = 1e-8,
    maxiter = 100,
    ncheck = 1,
    printout_every = typemax(Int),
    dtau_scaling = 1,
    halo = 1,
    viscosity_continuation::AbstractViscosityContinuation = NoViscosityContinuation(),
)
    grid.grid2d === grid.grid || throw(ArgumentError(
        "PseudoTransientSolver(::StaggeredGrid) requires a depth-averaged grid " *
        "(grid.grid2d === grid.grid, i.e. nz == 1); DIVA's vertical shear integral is " *
        "not yet ported (roadmaps/chmy.md, Phase 3). Build the grid without a " *
        "`layering` argument, e.g. `StaggeredGrid(T, lx, ly, dx, dy)`."))

    (; arch) = grid
    g = grid.grid2d
    T = eltype(g)
    acx() = _field(arch, g, NODE_ACX, T, halo)
    acy() = _field(arch, g, NODE_ACY, T, halo)
    return PseudoTransientSolver(
        T(ndim1), T(ndim2), T(ndim3), T(muB),
        T(theta_v), T(gamma), T(abstol), Int(maxiter), Int(ncheck),
        Int(printout_every), T(dtau_scaling),
        acx(), acy(), acx(), acy(), acx(), acy(), acx(), acy(),
        viscosity_continuation,
    )
end

"""
$(TYPEDSIGNATURES)

Solve the ice dynamics via a direct linear solver (e.g., sparse LU factorization).

# Improvements over `LegacyLinearMomentumSolver2D`:
 1. dynamics is a type parameter → dispatch on loop1!/loop2! is fully static, no need to thread a runtime `dynamics` argument through populate_vectors! layers.
 2. SparseMatrixCSC pre-allocated at construction (fixed sparsity pattern). The hot path writes directly to A.nzval via a precomputed COO→nzval index map, eliminating the sparse(Ai, Aj, Av) allocation on every solve.
 3. Single AI type parameter (i_idx and j_idx are always the same kind).
 4. VT/MT/PI type parameters for vector/matrix/perm arrays so the struct can hold GPU arrays (CuVector, CuSparseMatrix) without code changes. The populate_vectors! kernels are written with KernelAbstractions and run on whichever backend owns lsd.u.
"""
struct LinearMomentumSolver2D{
    DYN <: AbstractMomentumBalance,
    T   <: AbstractFloat,
    VT  <: AbstractVector,        # float vector type (u, u0, b)
    MT,                            # sparse matrix type (SparseMatrixCSC or CuSparseMatrix)
    PI  <: AbstractVector{Int},   # perm index vector type
    AI,
} <: AbstractMomentumSolver
    dynamics::DYN
    nx::Int
    ny::Int
    dxdx_::T
    dydy_::T
    dxdy_::T
    u::VT
    u0::VT
    b::VT
    A::MT
    perm::PI                      # COO fill order → A.nzval index
    i_idx::AI
    j_idx::AI
    solver_cache::Ref{Any}        # holds a cached direct solver (Nothing or backend-specific)
end

###############################################################
# Dispatch functions
###############################################################

# Functions to calculate velocity
function calc_F_integral(visc_eff, H_ice, f_ice, zeta_aa, n)
    # TODO: not yet implemented.
    error("calc_F_integral is not yet implemented")
end

# function velocity( solver::LinearSolver, dynamics::DIVA)
# end

function vertically_integrated_viscosity!(N, H, μ)
    N .= H .* μ
    return nothing
end

function stagger()
end

function fill_pattern!(Ai, Aj, nx, ny, i_idx, j_idx)
    k = 0
    for i in 1:nx, j in 1:ny
        im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
        nr = _ij2n_ux(i, j, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(ip1, j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(i,   j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(im1, j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(i,   jp1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(i,   jm1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(i,   j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(ip1, j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(ip1, jm1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(i,   jm1, nx, ny)
    end
    for i in 1:nx, j in 1:ny
        im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
        nr = _ij2n_uy(i, j, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(i,   jp1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(i,   j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(i,   jm1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(ip1, j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_uy(im1, j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(i,   jp1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(i,   j,   nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(im1, jp1, nx, ny)
        k+=1; Ai[k]=nr; Aj[k]=_ij2n_ux(im1, j,   nx, ny)
    end
    return nothing
end

function coo_to_nzval_idx(Ai, Aj, A::SparseMatrixCSC)
    perm = Vector{Int}(undef, length(Ai))
    for k in eachindex(Ai)
        c  = Aj[k]
        r  = Ai[k]
        lo = A.colptr[c]
        hi = A.colptr[c+1] - 1
        perm[k] = searchsortedfirst(A.rowval, r, lo, hi, Base.Order.Forward)
    end
    return perm
end

function LinearMomentumSolver2D(grid::RegularGrid, dynamics::DYN, backend = CPU()) where {DYN <: AbstractMomentumBalance}
    T      = eltype(grid.x)
    nx, ny = grid.nx, grid.ny
    dx, dy = T(grid.dx), T(grid.dy)
    dxdx_  = 1 / (dx * dx)
    dydy_  = 1 / (dy * dy)
    dxdy_  = 1 / (dx * dy)
    n_sprs = 18 * nx * ny  # 9 nonzeros/row × 2 equations (ux, uy)
    n_u    = 2 * nx * ny

    i_idx = PeriodicIndexing(1, nx)
    j_idx = PeriodicIndexing(1, ny)

    # Pattern and permutation are always computed on CPU (one-time cost).
    Ai_cpu   = zeros(Int, n_sprs)
    Aj_cpu   = zeros(Int, n_sprs)
    fill_pattern!(Ai_cpu, Aj_cpu, nx, ny, i_idx, j_idx)
    A_cpu    = sparse(Ai_cpu, Aj_cpu, ones(T, n_sprs), n_u, n_u)
    perm_cpu = coo_to_nzval_idx(Ai_cpu, Aj_cpu, A_cpu)

    # Allocate live arrays on the target backend.
    # For GPU (e.g. CUDABackend()), KernelAbstractions.zeros returns CuVector and
    # the sparse matrix should be adapted via CUDA.CUSPARSE.CuSparseMatrixCSC(A_cpu).
    u    = KernelAbstractions.zeros(backend, T,   n_u)
    u0   = KernelAbstractions.zeros(backend, T,   n_u)
    b    = KernelAbstractions.zeros(backend, T,   n_u)
    perm = KernelAbstractions.zeros(backend, Int, n_sprs)
    perm .= perm_cpu   # works for both CPU (no-op copy) and GPU (H→D transfer)

    # A stays as SparseMatrixCSC for the CPU default; for GPU, adapt before passing
    # to the inner constructor, e.g.:
    #   A = CUDA.CUSPARSE.CuSparseMatrixCSC(A_cpu)
    return LinearMomentumSolver2D(dynamics, nx, ny, dxdx_, dydy_, dxdy_, u, u0, b, A_cpu, perm, i_idx, j_idx, Ref{Any}(nothing))
end

@kernel function _assemble_ux!(nzval, perm, u0, b, nx, ny,
                               N, N_ab, ux, taud_acx, β_acx, β_acy,
                               dxdx_, dydy_, dxdy_, i_idx, j_idx, dynamics)
    i, j = @index(Global, NTuple)
    im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
    nr = _ij2n_ux(i, j, nx, ny)
    @inbounds u0[nr] = ux[i, j]
    @inbounds b[nr]  = taud_acx[i, j]
    v = loop1_coeffs(im1, i, ip1, jm1, j, jp1, dxdx_, dydy_, dxdy_,
                     N, N_ab, β_acx, β_acy, dynamics)
    k0 = 9 * ((i - 1) * ny + (j - 1))
    @inbounds for s in 1:9
        nzval[perm[k0 + s]] = v[s]
    end
end

@kernel function _assemble_uy!(nzval, perm, u0, b, nx, ny,
                               N, N_ab, uy, taud_acy, β_acx, β_acy,
                               dxdx_, dydy_, dxdy_, i_idx, j_idx, dynamics)
    i, j = @index(Global, NTuple)
    im1, ip1, jm1, jp1 = stencil(i, j, i_idx, j_idx)
    nr = _ij2n_uy(i, j, nx, ny)
    @inbounds u0[nr] = uy[i, j]
    @inbounds b[nr]  = taud_acy[i, j]
    v = loop2_coeffs(im1, i, ip1, jm1, j, jp1, dxdx_, dydy_, dxdy_,
                     N, N_ab, β_acx, β_acy, dynamics)
    k0 = 9 * nx * ny + 9 * ((i - 1) * ny + (j - 1))
    @inbounds for s in 1:9
        nzval[perm[k0 + s]] = v[s]
    end
end

"""
$(TYPEDSIGNATURES)

Assemble the linear system of the [`LinearMomentumSolver2D`](@ref) from the current dynamic
state `dyn_now`: fill the sparse matrix `A` (viscosity and basal-drag coefficients) and the
right-hand side `b` (driving stress), and seed the initial guess `u0` from the current
velocity. Mutates the solver buffers in place; the low-level method takes the unpacked fields
directly.
"""
function populate_vectors!(lsd::LinearMomentumSolver2D, dyn_now)
    (; N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy) = dyn_now
    populate_vectors!(lsd, N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy)
    return nothing
end

function populate_vectors!(
    lsd::LinearMomentumSolver2D,
    N, N_ab, ux, uy, taud_acx, taud_acy, β_acx, β_acy,
)
    (; A, perm, u0, b, i_idx, j_idx, nx, ny, dxdx_, dydy_, dxdy_, dynamics) = lsd
    backend  = get_backend(lsd.u)
    kern_ux  = _assemble_ux!(backend)
    kern_uy  = _assemble_uy!(backend)

    nzval = nonzeros(A)
    kern_ux(
        nzval, perm, u0, b, nx, ny,
        N, N_ab, ux, taud_acx, β_acx, β_acy,
        dxdx_, dydy_, dxdy_, i_idx, j_idx, dynamics;
        ndrange = (nx, ny),
    )
    kern_uy(
        nzval, perm, u0, b, nx, ny,
        N, N_ab, uy, taud_acy, β_acx, β_acy,
        dxdx_, dydy_, dxdy_, i_idx, j_idx, dynamics;
        ndrange = (nx, ny),
    )
    return nothing
end

function LinearSolve.LinearProblem(lsd::LinearMomentumSolver2D)
    return LinearProblem(lsd.A, lsd.b; u0 = lsd.u)
end

function velocity!(lsd::LinearMomentumSolver2D)
    if lsd.solver_cache[] === nothing
        F = lu(lsd.A)
        lsd.solver_cache[] = F
        ldiv!(lsd.u, F, lsd.b)
    else
        _velocity_cached!(lsd.solver_cache[], lsd)
    end
    return nothing
end

# Function barrier: typed on F so lu! and ldiv! dispatch statically.
function _velocity_cached!(F, lsd::LinearMomentumSolver2D)
    lu!(F, lsd.A)
    ldiv!(lsd.u, F, lsd.b)
    return nothing
end

function velocity!(ux, uy, lsd::LinearMomentumSolver2D)
    u = lsd.u
    (; nx, ny) = lsd
    @inbounds for i in 1:nx, j in 1:ny
        ux[i, j] = u[_ij2n_ux(i, j, nx, ny)]
        uy[i, j] = u[_ij2n_uy(i, j, nx, ny)]
    end
    return nothing
end
"""
$(TYPEDSIGNATURES)

Supertype for explicit time-integration methods.

# Subtypes
 - [`Euler`](@ref)
 - [`RungeKutta4`](@ref)
 - [`BogackiShampine32`](@ref)
 - [`Tsitouras54`](@ref)
 - [`RKL2`](@ref)
 - [`SSPRK33`](@ref)
 - [`SSPRK43`](@ref)
"""
abstract type AbstractIntegrationMethod end

# Fixed-step methods --------------------------------------------------------

"""
$(TYPEDSIGNATURES)

First-order, explicit (forward) Euler.

# Fields
 - `dt`: the fixed step size
 - `du`: scratch buffer for the derivative (same shape as the state)
"""
struct Euler{T, M} <: AbstractIntegrationMethod
    dt::T
    du::M
end
Euler(x::AbstractArray; dt) = Euler(eltype(x)(dt), similar(x))

"""
$(TYPEDSIGNATURES)

Explicit, 4th-order Runge–Kutta scheme.

# Fields
 - `dt`: the fixed step size
 - `du`: scratch buffer for the derivative (same shape as the state)
 - `k1..k4`: scratch buffers for the stage derivatives (same shape as the state)
"""
struct RungeKutta4{T, M} <: AbstractIntegrationMethod
    dt::T
    du::M
    k1::M
    k2::M
    k3::M
    k4::M
end
function RungeKutta4(x::AbstractArray; dt)
    RungeKutta4(eltype(x)(dt), similar(x), similar(x), similar(x), similar(x), similar(x))
end

# Adaptive methods ----------------------------------------------------------
#= NOTE: these are `mutable struct` so the error controller can write the next
# step size back into `dt` between calls. The scalar fields are isbits, so this
# costs nothing in type stability. The embedded error norm is computed by a fused
# reduction (see `_embedded_rms`), so no error-scratch buffer is stored.
=#

"""
$(TYPEDSIGNATURES)

Bogacki–Shampine 3(2) embedded pair. Third order with a second-order error
estimate.

# Fields
 - `dt`: the current step size
 - `du`: scratch buffer for the candidate state (same shape as the state)
 - `atol`, `rtol`: absolute and relative error tolerances
 - `dt_min`, `dt_max`: minimum and maximum allowed step sizes
 - `k1..k4`: stage derivatives

The weighted error norm is reduced in a single fused pass (see [`_embedded_rms`](@ref)),
so no error-scratch buffer is needed.
"""
mutable struct BogackiShampine32{T, M} <: AbstractIntegrationMethod
    dt::T
    du::M
    atol::T
    rtol::T
    dt_min::T
    dt_max::T
    k1::M
    k2::M
    k3::M
    k4::M
end
function BogackiShampine32(x::AbstractArray; dt, atol = 1e-6, rtol = 1e-3,
                           dt_min = 0, dt_max = Inf)
    T = eltype(x)
    BogackiShampine32(T(dt), similar(x), T(atol), T(rtol), T(dt_min), T(dt_max),
                      similar(x), similar(x), similar(x), similar(x))
end

"""
$(TYPEDSIGNATURES)

Tsitouras 5(4) embedded pair. Fifth order with a fourth-order error estimate.
Tsit5 is a **7-stage** FSAL method, so it needs seven stage buffers `k1..k7`
(your sketch had six); `du` is the candidate-state scratch. The weighted error norm
is reduced in a single fused pass (see [`_embedded_rms`](@ref)), so no error-scratch
buffer is needed.
"""
mutable struct Tsitouras54{T, M} <: AbstractIntegrationMethod
    dt::T
    du::M
    atol::T
    rtol::T
    dt_min::T
    dt_max::T
    k1::M
    k2::M
    k3::M
    k4::M
    k5::M
    k6::M
    k7::M
end
function Tsitouras54(x::AbstractArray; dt, atol = 1e-6, rtol = 1e-3,
                     dt_min = 0, dt_max = Inf)
    T = eltype(x)
    Tsitouras54(T(dt), similar(x), T(atol), T(rtol), T(dt_min), T(dt_max),
                similar(x), similar(x), similar(x), similar(x), similar(x),
                similar(x), similar(x))
end

# Stabilized explicit (super-time-stepping) ---------------------------------

"""
$(TYPEDSIGNATURES)

Second-order Runge–Kutta–Legendre super-time-stepping (RKL2; Meyer, Balsara &
Aslam, 2014). A stabilized *explicit* method for diffusion-dominated (stiff,
real-eigenvalue) problems such as the SIA thickness equation or vertical heat
conduction. It covers a superstep with `s` cheap stages whose stability region
stretches along the negative real axis (`Δt_stable ≈ Δt_expl · (s²+s−2)/4`), so it
takes far larger steps than an accuracy-controlled explicit pair (Tsit5/BS32) on
diffusive problems — while staying matrix-free (no Jacobian, no linear solve).

Unlike the embedded pairs, RKL2 is **fixed second order**; its adaptivity is in the
*stage count* `s`, chosen each superstep from the spectral radius. There is no embedded
accuracy estimate, so `dt` is a fixed target superstep, not error-adapted.

The spectral radius `ρ` can be supplied explicitly or estimated automatically each
superstep by a matrix-free power iteration (see [`estimate_spectral_radius!`](@ref)):

```julia
RKL2(x; dt = 1.0, spectral_radius = ρ)    # fixed ρ (cheapest if you know it analytically)
RKL2(x; dt = 1.0)                          # auto-estimate ρ each superstep
```

# Fields
 - `dt`: target superstep size.
 - `spectral_radius`: current `ρ` (upper bound on `|eigenvalues|` of the RHS Jacobian).
   In auto mode it is overwritten each superstep; otherwise it is the user-supplied value.
 - `auto`: whether to re-estimate `ρ` each superstep.
 - `safety`: safety factor applied to the estimated `ρ` (default 1.2).
 - `y0, f0, ym1, ym2, fm1`: stage buffers (same shape as the state).
 - `vest`: persistent eigenvector estimate that warm-starts the power iteration.
"""
mutable struct RKL2{T, M} <: AbstractIntegrationMethod
    dt::T
    spectral_radius::T
    auto::Bool
    safety::T
    y0::M
    f0::M
    ym1::M
    ym2::M
    fm1::M
    vest::M
end
function RKL2(x::AbstractArray; dt, spectral_radius = nothing, safety = 1.2)
    T    = eltype(x)
    auto = spectral_radius === nothing
    ρ0   = auto ? zero(T) : T(spectral_radius)
    RKL2(T(dt), ρ0, auto, T(safety),
         similar(x), similar(x), similar(x), similar(x), similar(x),
         fill!(similar(x), zero(T)))
end

# Strong-stability-preserving (SSP) methods ---------------------------------

"""
$(TYPEDSIGNATURES)

Three-stage, third-order strong-stability-preserving Runge–Kutta (SSPRK33, Shu–Osher
1988), SSP coefficient `C = 1`. Each stage is a convex combination of the previous stage
and a forward-Euler update, which guarantees that any monotonicity/positivity/TVD bound
satisfied by forward Euler is preserved — ideal for advective transport of a non-negative
quantity (ice thickness, tracers) with sharp fronts.

The admissible step is `C · dt_fe`, where `dt_fe` is the **forward-Euler-stable
(CFL-limited) step for the current state**, supplied externally (e.g.
`dt_fe = cfl · Δx / max|u|`). Update it whenever the state/velocities change.

# Fields
 - `dt_fe`: forward-Euler / CFL-limited stable step (user-set).
 - `us`: stage-state buffer (same shape as the state).
 - `L`: derivative buffer (same shape as the state).
"""
mutable struct SSPRK33{T, M} <: AbstractIntegrationMethod
    dt_fe::T
    us::M
    L::M
end
function SSPRK33(x::AbstractArray; dt_fe)
    T = eltype(x)
    SSPRK33(T(dt_fe), similar(x), similar(x))
end

"""
$(TYPEDSIGNATURES)

Four-stage, third-order strong-stability-preserving Runge–Kutta (SSPRK43; optimal
Spiteri–Ruuth 2002 pair), SSP coefficient `C = 2`. Same SSP guarantee as [`SSPRK33`](@ref),
but the extra stage doubles the SSP coefficient, so the work per unit stable time
(`C/stages = 0.5`) beats SSPRK33's `1/3` — usually the better choice.

Steps at `C · dt_fe = 2 · dt_fe`; see [`SSPRK33`](@ref) for the meaning of `dt_fe`.

# Fields
 - `dt_fe`: forward-Euler / CFL-limited stable step (user-set).
 - `us`: stage-state buffer (same shape as the state).
 - `L`: derivative buffer (same shape as the state).
"""
mutable struct SSPRK43{T, M} <: AbstractIntegrationMethod
    dt_fe::T
    us::M
    L::M
end
function SSPRK43(x::AbstractArray; dt_fe)
    T = eltype(x)
    SSPRK43(T(dt_fe), similar(x), similar(x))
end

# Integrator wrapper --------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Bundles a right-hand side `f`, the state `x`, parameters `p`, and an integration
`method`. The type parameters make every field concrete, so the in-place call
`f(du, x, p, t)` and the dispatch on `method` are statically resolved (no dynamic
dispatch, fully inlinable).

`f` must be an in-place, *non-capturing* function (a top-level `function` or a
callable struct) of the form `f(du, x, p, t)`.
"""
struct Integrator{F, M, P, I <: AbstractIntegrationMethod}
    f::F
    x::M
    p::P
    method::I
end

# Driver --------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Advance the integrator state in place over a macro-step `Δt_sync`, taking as many
internal sub-steps as the method requires. The right-hand side is treated as
autonomous over the window (boundary conditions are assumed frozen for the
duration of `Δt_sync`, consistent with the operator-splitting coupling in
`step!(::Simulation, …)`); stage times therefore run over a local `[0, Δt_sync]`.

Returns `nothing`; `integ.x` (and, for adaptive methods, `integ.method.dt`) are
mutated.
"""
# TODO: the macro-step assumes the RHS is autonomous over [0, Δt_sync] (BCs frozen).
#       If a component ever needs absolute time, thread a `t0` argument through here
#       and into the stage-time arguments of each `substep!`.
# TODO: optionally drive the step size from a CFL/stability limit (or
#       min(error-based, CFL-based)) rather than the embedded error estimate alone —
#       for ice the stability constraint usually binds before accuracy.
# NOTE: FSAL is intentionally not exploited (the state/params are mutated between
#       macro-steps by `couple!`, which would invalidate a carried-over first stage).
function step!(integ::Integrator, Δt_sync)
    (; f, x, p, method) = integ
    t  = zero(Δt_sync)
    while t < Δt_sync
        t = substep!(f, x, p, method, t, Δt_sync)
    end
    return nothing
end

# Fixed-step sub-steps ------------------------------------------------------

function substep!(f::F, x, p, m::Euler, t, t_end) where {F}
    h = min(m.dt, t_end - t)
    f(m.du, x, p, t)
    @. x += h * m.du
    return t + h
end

function substep!(f::F, x, p, m::RungeKutta4, t, t_end) where {F}
    h = min(m.dt, t_end - t)
    f(m.k1, x, p, t)
    @. m.du = x + (h / 2) * m.k1;  f(m.k2, m.du, p, t + h / 2)
    @. m.du = x + (h / 2) * m.k2;  f(m.k3, m.du, p, t + h / 2)
    @. m.du = x + h * m.k3;        f(m.k4, m.du, p, t + h)
    @. x += (h / 6) * (m.k1 + 2 * m.k2 + 2 * m.k3 + m.k4)
    return t + h
end

# Adaptive sub-steps --------------------------------------------------------
#
# TODO (applies to both BogackiShampine32 and Tsitouras54 substeps below):
#  - Final-substep dt leak: when the last substep is clamped to land exactly on
#    `t_end`, that smaller `h` is written back into `m.dt` and carries into the next
#    macro-step's first step. Benign (self-corrects in a step or two) but could be
#    fixed by not updating `m.dt` when the step was limited by `t_end` rather than error.
#  - Rejection-loop guard: with `dt_min == 0` a pathological RHS could loop forever.
#    Add a max-attempts cap (and/or error out) once `h` underflows `dt_min`.
#  - Replace the elementary I-controller (`fac = safety·err^expo`) with a PI controller
#    for smoother step-size sequences on stiff/non-smooth problems.

"""
$(TYPEDSIGNATURES)

RMS of the embedded error, `sqrt(Σ wᵢ² / n)`, with element-wise weighted error
`wᵢ = h·Σⱼ eⱼ·kⱼ[i] / (atol + rtol·max(|x[i]|, |du[i]|))`. `sum(abs2, ·)` reduces a *lazy*
broadcast in a single fused pass — one GPU kernel, no materialized temporary — so the
adaptive methods carry no error buffer.

`h`, the tolerances, and the error weights are passed as arguments (not captured) so the
inner broadcast closure holds no boxed loop variable; capturing the reassigned step size
`h` directly would box it and allocate on every step. Two arities serve the 4-stage
(BS32) and 7-stage (Tsit5) pairs.
"""
@inline function _embedded_rms(h, x, du, atol, rtol,
                               k1, k2, k3, k4, e1, e2, e3, e4)
    werr = Broadcast.instantiate(Broadcast.broadcasted(
        x, du, k1, k2, k3, k4) do xi, dui, s1, s2, s3, s4
        h * (e1 * s1 + e2 * s2 + e3 * s3 + e4 * s4) / (atol + rtol * max(abs(xi), abs(dui)))
    end)
    return sqrt(sum(abs2, werr) / length(x))
end
@inline function _embedded_rms(h, x, du, atol, rtol,
                               k1, k2, k3, k4, k5, k6, k7,
                               e1, e2, e3, e4, e5, e6, e7)
    werr = Broadcast.instantiate(Broadcast.broadcasted(
        x, du, k1, k2, k3, k4, k5, k6, k7) do xi, dui, s1, s2, s3, s4, s5, s6, s7
        h * (e1 * s1 + e2 * s2 + e3 * s3 + e4 * s4 + e5 * s5 + e6 * s6 + e7 * s7) /
            (atol + rtol * max(abs(xi), abs(dui)))
    end)
    return sqrt(sum(abs2, werr) / length(x))
end

"""
$(TYPEDSIGNATURES)

Butcher tableau for the method, returned as a `NamedTuple` of coefficients at the
method's precision `T` plus `expo`, the step-controller exponent `-1/(p̂+1)` where
`p̂` is the embedded order. Dispatching on the method *type* keeps these compile-time
constants out of the (mutable) instance struct, so they are inlined and
constant-folded into `substep!` rather than loaded from memory.
"""
@inline function tableau(::BogackiShampine32{T}) where {T}
    return (
        a21 = T(1//2), a32 = T(3//4),                                  # stage / node coeffs
        b1  = T(2//9),  b2 = T(1//3),  b3 = T(4//9),                   # 3rd-order weights
        e1  = T(-5//72), e2 = T(1//12), e3 = T(1//9), e4 = T(-1//8),   # error weights (bᵢ − b̂ᵢ)
        expo = -one(T) / 3,                                            # embedded order 2
    )
end

@inline function tableau(::Tsitouras54{T}) where {T}
    return (
        # Tsit5 tableau (Tsitouras, 2011)
        c2 = T(0.161), c3 = T(0.327), c4 = T(0.9), c5 = T(0.9800255409045097),
        a21 = T(0.161),
        a31 = T(-0.008480655492356989), a32 = T(0.335480655492357),
        a41 = T(2.8971530571054935),    a42 = T(-6.359448489975075),  a43 = T(4.3622954328695815),
        a51 = T(5.325864828439257),     a52 = T(-11.748883564062828), a53 = T(7.4955393428898365), a54 = T(-0.09249506636175525),
        a61 = T(5.86145544294642),      a62 = T(-12.92096931784711),  a63 = T(8.159367898576159),  a64 = T(-0.071584973281401), a65 = T(-0.028269050394068383),
        # 5th-order weights bᵢ == a7ᵢ
        b1 = T(0.09646076681806523), b2 = T(0.01), b3 = T(0.4798896504144996),
        b4 = T(1.379008574103742),   b5 = T(-3.290069515436081), b6 = T(2.324710524099774),
        # error weights bᵢ − b̂ᵢ
        bt1 = T(-0.001780011052226), bt2 = T(-0.000816434459657), bt3 = T(0.007880878010262),
        bt4 = T(-0.144711007173263), bt5 = T(0.582357165452555),  bt6 = T(-0.458082105929187), bt7 = T(0.015151515151515),
        expo = -one(T) / 5,                                          # embedded order 4
    )
end

function substep!(f::F, x, p, m::BogackiShampine32, t, t_end) where {F}
    T = eltype(x)
    (; a21, a32, b1, b2, b3, e1, e2, e3, e4, expo) = tableau(m)
    safety = T(0.9); minfac = T(0.2); maxfac = T(10)
    atol, rtol = m.atol, m.rtol
    remaining  = t_end - t
    h = min(m.dt, remaining)
    while true
        f(m.k1, x, p, t)
        @. m.du = x + h * a21 * m.k1;  f(m.k2, m.du, p, t + a21 * h)
        @. m.du = x + h * a32 * m.k2;  f(m.k3, m.du, p, t + a32 * h)
        @. m.du = x + h * (b1 * m.k1 + b2 * m.k2 + b3 * m.k3)
        f(m.k4, m.du, p, t + h)              # candidate at du is the 3rd-order solution
        err = _embedded_rms(h, x, m.du, atol, rtol, m.k1, m.k2, m.k3, m.k4, e1, e2, e3, e4)
        fac = clamp(safety * err^expo, minfac, maxfac)
        if err <= one(T) || h <= m.dt_min
            @. x = m.du
            m.dt = clamp(h * fac, m.dt_min, m.dt_max)
            return t + h
        end
        h = clamp(h * fac, m.dt_min, remaining)
    end
end

function substep!(f::F, x, p, m::Tsitouras54, t, t_end) where {F}
    T = eltype(x)
    (; c2, c3, c4, c5,
       a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65,
       b1, b2, b3, b4, b5, b6,
       bt1, bt2, bt3, bt4, bt5, bt6, bt7, expo) = tableau(m)
    safety = T(0.9); minfac = T(0.2); maxfac = T(10)
    atol, rtol = m.atol, m.rtol

    remaining = t_end - t
    h = min(m.dt, remaining)
    while true
        f(m.k1, x, p, t)
        @. m.du = x + h * (a21 * m.k1);                                              f(m.k2, m.du, p, t + c2 * h)
        @. m.du = x + h * (a31 * m.k1 + a32 * m.k2);                                 f(m.k3, m.du, p, t + c3 * h)
        @. m.du = x + h * (a41 * m.k1 + a42 * m.k2 + a43 * m.k3);                    f(m.k4, m.du, p, t + c4 * h)
        @. m.du = x + h * (a51 * m.k1 + a52 * m.k2 + a53 * m.k3 + a54 * m.k4);       f(m.k5, m.du, p, t + c5 * h)
        @. m.du = x + h * (a61 * m.k1 + a62 * m.k2 + a63 * m.k3 + a64 * m.k4 + a65 * m.k5)
        f(m.k6, m.du, p, t + h)
        # 5th-order candidate solution
        @. m.du = x + h * (b1 * m.k1 + b2 * m.k2 + b3 * m.k3 + b4 * m.k4 + b5 * m.k5 + b6 * m.k6)
        f(m.k7, m.du, p, t + h)
        err = _embedded_rms(h, x, m.du, atol, rtol,
                            m.k1, m.k2, m.k3, m.k4, m.k5, m.k6, m.k7,
                            bt1, bt2, bt3, bt4, bt5, bt6, bt7)
        fac = clamp(safety * err^expo, minfac, maxfac)
        if err <= one(T) || h <= m.dt_min
            @. x = m.du
            m.dt = clamp(h * fac, m.dt_min, m.dt_max)
            return t + h
        end
        h = clamp(h * fac, m.dt_min, remaining)
    end
end

# RKL2 Legendre weights b_j (b_0 = b_1 = 1/3; the j≥2 formula also gives b_2 = 1/3).
@inline _rkl2_b(j, ::Type{T}) where {T} =
    j ≥ 2 ? (T(j)^2 + T(j) - 2) / (2 * T(j) * (T(j) + 1)) : T(1//3)

"""
$(TYPEDSIGNATURES)

Matrix-free estimate of the spectral radius `ρ = max|λ(∂f/∂x)|` at the current state,
via a nonlinear power method (the scheme used by RKC/ROCK). Jacobian–vector products
are approximated by finite differences `J·v ≈ (f(x+εv) − f(x))/ε`, so it needs only RHS
evaluations and works with any `f`, on CPU or GPU.

Assumes the caller has already set `m.y0 == x` and `m.f0 == f(x, p, t)`. Writes
`m.safety · ρ` into `m.spectral_radius` and stores the dominant eigenvector estimate in
`m.vest` to warm-start (and thus shorten) the next call. Returns the estimate. Buffers
`m.ym1`/`m.fm1` are used as scratch and overwritten.
"""
# TODO: clustered spectra (e.g. a real diffusion operator, whose top eigenvalues differ
#       by <1%) converge slowly here, so the *cold-start* estimate can under-resolve —
#       and underestimating ρ is the unsafe direction for stability. Mitigations to add:
#        - prefer an analytic bound when the operator provides one (e.g. SIA diffusion:
#          ρ ≈ 2D(1/Δx² + 1/Δy²)); fall back to this power method only otherwise;
#        - raise `maxiter` / tighten the stopping tolerance for the first (cold) call;
#        - inflate `safety` for cold starts.
# TODO: complex-dominated spectra (strong advection) violate RKL2's real-axis stability
#       assumption; the magnitude returned here is then not a sufficient stability bound.
#       Use an SSP/IMEX integrator for advection-dominated components instead.
function estimate_spectral_radius!(f::F, p, t, m::RKL2; maxiter = 50) where {F}
    T = eltype(m.y0)
    yn, fn, v, fv, vest = m.y0, m.f0, m.ym1, m.fm1, m.vest
    uround = eps(T)

    nrm_y = norm(yn)
    # Cold start from a *broadband* random direction so the iteration can see the
    # dominant (often high-frequency) eigenvector; a smooth/modal guess like f(x) can
    # stay trapped in a sub-dominant eigenspace. Warm starts reuse the converged vest.
    norm(vest) == 0 && randn!(vest)
    nrm_v = norm(vest)

    dynrm = nrm_y != 0 ? nrm_y * sqrt(uround) : uround
    @. v = yn + vest * (dynrm / nrm_v)              # perturb x by magnitude dynrm along vest

    σ = zero(T)
    for iter in 1:maxiter
        f(fv, v, p, t)
        @. fv = fv - fn                              # ≈ J·(v − yn)
        dfnrm = norm(fv)
        σ_old = σ
        σ = dfnrm / dynrm
        (iter ≥ 2 && abs(σ - σ_old) ≤ T(0.01) * max(σ, σ_old)) && break
        dfnrm == 0 && break
        @. v = yn + fv * (dynrm / dfnrm)            # next power-iteration direction
        copyto!(vest, fv)                            # remember it for the warm start
    end
    m.spectral_radius = m.safety * σ
    return m.spectral_radius
end

function substep!(f::F, x, p, m::RKL2, t, t_end) where {F}
    T = eltype(x)
    h = min(m.dt, t_end - t)

    # Y_0 and F_0 = M(Y_0)
    copyto!(m.y0, x)
    f(m.f0, m.y0, p, t)
    # TODO: re-estimating ρ every superstep costs (a few) extra f-evals; for cheap
    #       warm-started estimates that's fine, but consider re-estimating only every
    #       N supersteps (ρ varies slowly) and reusing m.spectral_radius in between.
    m.auto && estimate_spectral_radius!(f, p, t, m) # refresh ρ from the current state

    # Stage count for stability: h ≤ Δt_expl·(s²+s−2)/4 with Δt_expl = 2/ρ.
    # Invert for the smallest integer s (≥ 2).
    # TODO: `s` is uncapped — for extreme stiffness (h·ρ ≫ 1) it grows like √(hρ) and can
    #       become very large. Cap `s` at some s_max and take multiple supersteps instead.
    ρ = m.spectral_radius
    r = ρ > 0 ? h * ρ / 2 : zero(T)                 # = h / Δt_expl
    s = max(2, ceil(Int, (sqrt(T(9) + 16 * r) - 1) / 2))
    w1 = 4 / (T(s)^2 + T(s) - 2)

    μ̃1 = _rkl2_b(1, T) * w1
    @. m.ym1 = m.y0 + μ̃1 * h * m.f0                 # Y_1
    copyto!(m.ym2, m.y0)                             # Y_0 (= Y_{j-2} at j = 2)

    for j in 2:s
        bj   = _rkl2_b(j,     T)
        bjm1 = _rkl2_b(j - 1, T)
        bjm2 = _rkl2_b(j - 2, T)
        μj  = T(2j - 1) / T(j) * bj / bjm1
        νj  = -T(j - 1) / T(j) * bj / bjm2
        μ̃j  = μj * w1
        γ̃j  = -(1 - bjm1) * μ̃j                      # -a_{j-1} μ̃_j,  a_{j-1} = 1 - b_{j-1}
        cj0 = 1 - μj - νj
        f(m.fm1, m.ym1, p, t)                       # M(Y_{j-1}); f autonomous over the window
        @. m.ym2 = μj * m.ym1 + νj * m.ym2 + cj0 * m.y0 + μ̃j * h * m.fm1 + γ̃j * h * m.f0
        m.ym1, m.ym2 = m.ym2, m.ym1                 # rotate: ym1 = Y_j, ym2 = Y_{j-1} (O(1) swap)
    end
    copyto!(x, m.ym1)                               # Y_s
    return t + h
end

# SSP sub-steps -------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

SSP coefficient `C`: the method preserves forward Euler's monotonicity/positivity/TVD
bound for steps up to `C · dt_fe`. `C = 1` for [`SSPRK33`](@ref), `C = 2` for [`SSPRK43`](@ref).
"""
ssp_coefficient(::SSPRK33{T}) where {T} = one(T)
ssp_coefficient(::SSPRK43{T}) where {T} = T(2)

# TODO: SSPRK43 admits an embedded 2nd-order companion for error-based adaptivity; could
#       add it (mirroring the BS32/Tsit5 controller) if accuracy control is ever wanted
#       on top of the CFL/SSP step limit.

function substep!(f::F, x, p, m::SSPRK33, t, t_end) where {F}
    T = eltype(x)
    h = min(ssp_coefficient(m) * m.dt_fe, t_end - t)
    # Shu–Osher convex-combination form (SSP-preserving for h ≤ C·dt_fe).
    f(m.L, x, p, t)
    @. m.us = x + h * m.L                                       # u⁽¹⁾
    f(m.L, m.us, p, t)
    @. m.us = T(3//4) * x + T(1//4) * m.us + T(1//4) * h * m.L  # u⁽²⁾
    f(m.L, m.us, p, t)
    @. x = T(1//3) * x + T(2//3) * m.us + T(2//3) * h * m.L     # uⁿ⁺¹
    return t + h
end

function substep!(f::F, x, p, m::SSPRK43, t, t_end) where {F}
    T = eltype(x)
    h = min(ssp_coefficient(m) * m.dt_fe, t_end - t)
    f(m.L, x, p, t)
    @. m.us = x + T(1//2) * h * m.L                             # u⁽¹⁾
    f(m.L, m.us, p, t)
    @. m.us = m.us + T(1//2) * h * m.L                          # u⁽²⁾
    f(m.L, m.us, p, t)
    @. m.us = T(2//3) * x + T(1//3) * m.us + T(1//6) * h * m.L  # u⁽³⁾
    f(m.L, m.us, p, t)
    @. x = m.us + T(1//2) * h * m.L                             # uⁿ⁺¹
    return t + h
end

# Example usage:
#
#   method = Tsitouras54(m0; dt = 1.0, atol = 1e-6, rtol = 1e-3)
#   integ  = Integrator(dmdt!, m0, p, method)   # dmdt!(du, x, p, t)
#   step!(integ, Δt_sync)
#
# Verify type stability with:  using Test; @inferred step!(integ, Δt_sync)

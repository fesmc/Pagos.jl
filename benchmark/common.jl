# ---------------------------------------------------------------------------
# Shared fixtures for the tracked benchmark suite.
#
# Everything a timed function reads is built here, once, from closed-form expressions —
# never `rand`. Two commits' numbers are only comparable if the inputs are bit-identical,
# and several kernels on the hot path branch on their data (`node_active` on the mask, the
# viscosity regularisation floor, the `H > 0` guards in `dotvel!`), so random inputs would
# let the *fraction of cells taking each branch* drift between runs.
#
# The geometry is the uniform ice slab of `test/mechanics/pseudotransient_staggered.jl`
# (H0 = 1000 m, μ0 = 1e5, β0 = 1e4, α = 1e-3): uniform fields make the membrane-stress
# divergence vanish, so the PT solve has a closed-form answer and converges in ~10
# iterations, deterministically. That determinism is the point — a converged-solve
# benchmark whose iteration count wanders is not a measurement of anything. Kernel cost is
# independent of the *values* here (same node classes, same sweep, same branch taken
# everywhere), so the uniformity costs nothing on the per-kernel entries.
# ---------------------------------------------------------------------------

using Pagos
using BenchmarkTools
using KernelAbstractions

# Grid the tracked suite reports on. One size: this is a regression tracker, not a scaling
# study — resolution sweeps belong in `basics/`.
const BENCH_T  = Float64
const BENCH_NX = 256
const BENCH_NY = 256
const BENCH_NZ = 11
const BENCH_DX = 5.0e3     # m

# Slab parameters, matching the analytic case the tests pin.
const BENCH_H0    = 1.0e3     # m
const BENCH_MU0   = 1.0e5     # Pa yr
const BENCH_BETA0 = 1.0e4     # Pa yr m⁻¹
const BENCH_ALPHA = 1.0e-3    # surface slope
const BENCH_A0    = 1.0e-16   # rate factor, Pa⁻ⁿ yr⁻¹

"""
    gpu_requested()

Whether the gated GPU group should be built: `PAGOS_BENCH_GPU=1` *and* a functional CUDA
device. Off by default so a CPU baseline stays comparable across machines, and so the
suite runs at all on a box without a GPU.
"""
function gpu_requested()
    get(ENV, "PAGOS_BENCH_GPU", "0") in ("1", "true", "yes") || return false
    return Base.find_package("CUDA") !== nothing
end

"""
    sync!(backend)

Wait for every kernel launched on `backend` to finish. **Must** close every timed region:
Pagos' kernels launch asynchronously on GPU, so without this a GPU benchmark measures
launch overhead rather than execution. A no-op on the CPU backend, where kernels are
already synchronous, so the same benchmarkable body serves both.
"""
@inline sync!(backend) = KernelAbstractions.synchronize(backend)

"""
    to_backend(backend, a)

Copy the host array `a` onto `backend`, returning an array of that backend's own type.
"""
function to_backend(backend, a::AbstractArray)
    backend isa KernelAbstractions.CPU && return copy(a)
    d = KernelAbstractions.allocate(backend, eltype(a), size(a)...)
    copyto!(d, a)
    return d
end

"""
    fill_field!(f, grid, fun)

Fill `f`'s interior *and* both halo rings from `fun(x, y, ζ)` at `f`'s own node class.
The halo matters: `Chmy.set!`/`setdata!` touch the interior only, which leaves every
boundary stencil reading an unset ghost cell — harmless for a correctness-free timing run,
except that unset memory can hold `NaN`, and `NaN` changes the branch a mask or a `max`
guard takes. Same helper the Chmy-native tests use (`test/test_helpers/chmy.jl`).
"""
function fill_field!(f, grid, fun)
    loc = location(f)
    ni, nj, nk = size(interior(f))
    for k in -1:(nk + 2), j in -1:(nj + 2), i in -1:(ni + 2)
        x, y, ζ = coord(grid, loc, i, j, k)
        f[i, j, k] = fun(x, y, ζ)
    end
    return f
end

"""
    bench_fixture(; backend, nx, ny, nz, dx, T)

Build the shared state every `dynamics/` benchmark runs against: the grid and its
[`Runtime`](@ref), a `MechanicState`/`MaterialState` filled with the uniform slab, a
`MaterialState` for the 3D deviatoric stress, and the `IceMask` a production solve uses.

Returned as a `NamedTuple` so a benchmark file destructures only what it needs.

!!! note "The mask is `IceMask`, not `NoMask`"
    Every masked kernel carries a `node_active` branch that a real solve pays for; timing
    the `NoMask` path would track code the solver does not run. The slab is fully
    ice-covered, so the branch is uniformly taken — cost, not divergence, is what this
    measures.
"""
function bench_fixture(; backend = KernelAbstractions.CPU(),
                       nx = BENCH_NX, ny = BENCH_NY, nz = BENCH_NZ,
                       dx = BENCH_DX, T = BENCH_T)
    layering = CorrectedVerticalLayering(T, QuadraticSigmaTransform(T, nz))
    grid = StaggeredGrid(backend, T, T(nx * dx), T(ny * dx), T(dx), T(dx), layering)
    rt   = Runtime(grid)
    cst  = Constants{T}()

    mech = MechanicState(grid)
    mat  = MaterialState(grid)

    H0, μ0, β0, α, A0 = T(BENCH_H0), T(BENCH_MU0), T(BENCH_BETA0), T(BENCH_ALPHA), T(BENCH_A0)

    # Geometry and material, uniform in the horizontal; the surface carries the slope that
    # drives the flow. Viscosity varies with ζ so the DIVA integrals see a real column
    # (a depth-uniform µ makes F₁ exact by construction — see the viscosity-integral tests
    # — which would flatter the quadrature without changing its cost).
    fill_field!(mech.topography.thickness, rt.grid2d, (x, y, ζ) -> H0)
    fill_field!(mech.topography.surface,   rt.grid2d, (x, y, ζ) -> H0 - α * x)
    fill_field!(mech.material.viscosity_depthaveraged, rt.grid2d, (x, y, ζ) -> μ0)
    fill_field!(mech.material.viscosity, rt.grid, (x, y, ζ) -> μ0 * (1 + ζ))
    fill_field!(mech.material.rate_factor_depthaveraged, rt.grid2d, (x, y, ζ) -> A0)
    fill_field!(mech.material.rate_factor, rt.grid, (x, y, ζ) -> A0)
    fill_field!(mech.friction.beta,     rt.grid2d, (x, y, ζ) -> β0)
    fill_field!(mech.friction.beta_eff, rt.grid2d, (x, y, ζ) -> β0)
    fill_field!(mat.eta_ice, rt.grid, (x, y, ζ) -> μ0 * (1 + ζ))

    # A non-zero, non-uniform velocity: the strain-rate and stress kernels are timed on
    # their own, outside any solve, and a zero velocity field would send the effective
    # strain rate straight to its regularisation floor — a branch the solver only takes on
    # iteration 1.
    u0 = T(50)
    fill_field!(mech.velocity.depthaverage_x, rt.grid2d, (x, y, ζ) -> u0 * (1 + 1e-6x))
    fill_field!(mech.velocity.depthaverage_y, rt.grid2d, (x, y, ζ) -> u0 * 1e-6y)
    fill_field!(mech.velocity.base_x, rt.grid2d, (x, y, ζ) -> u0 / 2)
    fill_field!(mech.velocity.base_y, rt.grid2d, (x, y, ζ) -> zero(T))
    fill_field!(mech.velocity.x, rt.grid, (x, y, ζ) -> u0 * (1 + ζ) / 2)
    fill_field!(mech.velocity.y, rt.grid, (x, y, ζ) -> zero(T))

    # The mask a production solve carries, built the way `docs/src/examples/ais-momentum.jl`
    # builds it: from a `TopographicState`, not hand-assembled.
    topo = TopographicState(grid)
    fill_field!(topo.thickness.ice, rt.grid2d, (x, y, ζ) -> H0)
    fill_field!(topo.mask.is_grounded, rt.grid2d, (x, y, ζ) -> true)
    icemasks!(topo, rt)
    momentum_mask!(topo, rt)
    mask = IceMask(topo.mask.is_momentum_solved)

    snapshot = map(f -> copy(asarray(f)), _momentum_inputs(mech))

    return (; grid, rt, cst, mech, mat, mask, topo, snapshot, backend, T, nx, ny, nz, dx)
end

# The fields the momentum path *reads* and some other benchmark in the suite *writes*.
# `beta_eff_diva!` rewrites β_eff (by a factor of ~30 on this slab), `update_viscosity!`
# and `depthaverage!` rewrite the viscosities, `pseudo_transient!` rewrites the velocity —
# all inputs to the solves. Benchmarks within a `BenchmarkGroup` run in dictionary order,
# which is not the source order and is not stable across Julia versions, so without an
# explicit restore a solve's iteration count would depend on what happened to run before
# it. That is precisely the kind of hidden state that makes two commits' numbers
# incomparable.
_momentum_inputs(mech) = (
    mech.velocity.depthaverage_x, mech.velocity.depthaverage_y,
    mech.velocity.base_x, mech.velocity.base_y,
    mech.friction.beta_eff,
    mech.material.viscosity, mech.material.viscosity_depthaveraged,
    mech.material.viscosity_integral_1, mech.material.viscosity_integral_2,
)

"""
    restore_momentum!(fx)

Put every momentum input back to what `bench_fixture` built, from the snapshot taken there.
Belongs in a benchmark's `setup`, never in its timed body.

A `copyto!` per field (halos included), not a re-fill: restoring must be cheap enough to
run before every sample, and it must reproduce the *bytes* the fixture started with rather
than re-deriving them.
"""
function restore_momentum!(fx)
    for (f, saved) in zip(_momentum_inputs(fx.mech), fx.snapshot)
        copyto!(asarray(f), saved)
    end
    return nothing
end

"""
    reset_velocity!(fx)

Zero the depth-averaged velocity the PT solver iterates on — the initial guess every solve
benchmark starts from. Without it, sample *n* would start from sample *n-1*'s converged
answer and time a solve that had already happened.
"""
function reset_velocity!(fx)
    setdata!(fx.mech.velocity.depthaverage_x, zero(fx.T))
    setdata!(fx.mech.velocity.depthaverage_y, zero(fx.T))
    return nothing
end

# ---------------------------------------------------------------------------
# Shared fixture and helpers for the GPU probes in this directory.
#
# The geometry is the same closed-form ice slab the tracked suite uses
# (`benchmark/common.jl`), at a resolution large enough for a GPU to be bandwidth-bound
# rather than launch-bound. No data file, no `rand`: two runs are comparable, and the
# kernels' cost is independent of the values anyway (uniform slab ⇒ every `node_active`
# branch taken the same way at every node, which is what a fully ice-covered interior of a
# real ice sheet looks like too).
#
# `nz = 2` is the smallest column DIVA accepts (`_check_momentum_grid` rejects `nz == 1`),
# and under the default `NoDIVUpdate` the PT loop touches only depth-integrated fields —
# so the loop cost is independent of `nz` and the small column keeps a Float64 state inside
# 8 GB.
# ---------------------------------------------------------------------------

using Pagos
using Pagos.KernelAbstractions
using Pagos.Chmy
using Pagos.Chmy.Architectures: get_backend, heuristic_groupsize
using Pagos: NODE_AA, NODE_AB, NODE_ACX, NODE_ACY
using CUDA
using Printf

const GPU_NX = 760
const GPU_NY = 760
const GPU_NZ = 2
const GPU_DX = 8.0e3      # m, the AIS restarts' 8 km spacing

const SLAB_H0    = 1.0e3     # m
const SLAB_MU0   = 1.0e5     # Pa yr
const SLAB_BETA0 = 1.0e4     # Pa yr m⁻¹
const SLAB_ALPHA = 1.0e-3    # surface slope

# Solver configuration these probes measure — the one `docs/src/examples/ais-momentum/`
# runs. `abstol` and `maxiter` are overridden per call so a probe times a *fixed* iteration
# count rather than a convergence path.
const GPU_SOLVER_KWARGS = (
    ncheck = 20, printout_every = typemax(Int),
    pseudo_timestep = GershgorinPseudoTimeStep(cfl = 0.8),
    convergence = ScaledResidual(),
    friction_update = ActiveFrictionUpdate(),
    tuning = AutotunedDynamicRelaxation(cadence = 50),
)

"""
    fill_field!(f, grid, fun)

Fill `f`'s interior *and* both halo rings from `fun(x, y, ζ)` at `f`'s own node class.
Element-by-element, so it runs on a **CPU** field only — see [`gpu_fixture`](@ref).
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
    gpu_fixture(; T, nx, ny, nz, dx)

Build the slab `MechanicState`, mask and `Runtime` **on the GPU**.

Built on the CPU first and moved across in one `Adapt.adapt(CuArray, ...)`: a `CuArray`
disallows scalar indexing, so [`fill_field!`](@ref)'s element-by-element writes cannot run
on a device field. This is the same two-step the `ais-momentum` example uses, and the
reason the tracked suite's `bench_fixture(; backend = CUDABackend())` cannot be reused
here.
"""
function gpu_fixture(; T = Float64, nx = GPU_NX, ny = GPU_NY, nz = GPU_NZ, dx = GPU_DX)
    arch = Arch(CUDABackend())
    layering() = CorrectedVerticalLayering(T, QuadraticSigmaTransform(T, nz))
    mk(a) = StaggeredGrid(a..., T, T(nx * dx), T(ny * dx), T(dx), T(dx), layering())
    grid_gpu = mk((arch,))
    grid_cpu = mk(())
    rt_cpu = Runtime(grid_cpu)

    H0, μ0, β0, α = T(SLAB_H0), T(SLAB_MU0), T(SLAB_BETA0), T(SLAB_ALPHA)
    mech = MechanicState(grid_cpu)
    fill_field!(mech.topography.thickness, rt_cpu.grid2d, (x, y, ζ) -> H0)
    fill_field!(mech.topography.surface,   rt_cpu.grid2d, (x, y, ζ) -> H0 - α * x)
    fill_field!(mech.material.viscosity_depthaveraged, rt_cpu.grid2d, (x, y, ζ) -> μ0)
    fill_field!(mech.material.viscosity, rt_cpu.grid, (x, y, ζ) -> μ0 * (1 + ζ))
    fill_field!(mech.friction.beta,     rt_cpu.grid2d, (x, y, ζ) -> β0)
    fill_field!(mech.friction.beta_eff, rt_cpu.grid2d, (x, y, ζ) -> β0)
    setdata!(mech.velocity.depthaverage_x, zero(T))
    setdata!(mech.velocity.depthaverage_y, zero(T))

    topo = TopographicState(grid_cpu)
    fill_field!(topo.thickness.ice, rt_cpu.grid2d, (x, y, ζ) -> H0)
    fill_field!(topo.mask.is_grounded, rt_cpu.grid2d, (x, y, ζ) -> true)
    icemasks!(topo, rt_cpu)
    momentum_mask!(topo, rt_cpu)

    topo_gpu = Pagos.Adapt.adapt(CuArray, topo)
    fx = (; grid = grid_gpu, rt = Runtime(grid_gpu), cst = Constants{T}(),
          mech = Pagos.Adapt.adapt(CuArray, mech),
          mask = IceMask(topo_gpu.mask.is_momentum_solved), T, nx, ny, nz)
    GC.gc(); CUDA.reclaim()
    return fx
end

"""
    bench(f, label; n = 50, reps = 5)

Time `f` and report µs per call: the **minimum** over `reps` batches of `n` calls. GPU
clocks drift under sustained load (badly so on a Max-Q part), and the minimum is the figure
that reproduces across runs; a mean tracks the thermal state of the machine instead.
"""
function bench(f, label = nothing; n = 50, reps = 5)
    for _ in 1:3; f(); end
    CUDA.synchronize()
    us = minimum(_ -> 1e6 * CUDA.@elapsed(begin
                                              for _ in 1:n; f(); end
                                          end) / n, 1:reps)
    label === nothing || @printf("  %-46s %8.1f us\n", label, us)
    return us
end

# ---------------------------------------------------------------------------
# FlatLauncher — see `layout_and_launch.jl` for what it is worth and why.
# ---------------------------------------------------------------------------

"""
    FlatLauncher(arch, grid2d; sync = true)

A `Chmy.Launcher` replacement for **depth-integrated** grids, differing in two ways:

 1. worksize `(nx + 2, ny + 2, 1)` with `Offset(-1, -1, 0)`, against Chmy's
    `size(grid, Center()) .+ 2` = `(nx + 2, ny + 2, 3)` with `Offset(-1)`. `grid2d`'s z axis
    has extent 1, so the two extra planes Chmy's rule adds are a halo ring in a dimension
    that has none — every depth-integrated kernel does 3× the work it needs to.
 2. `sync = false` drops the `KernelAbstractions.synchronize(backend)` Chmy runs after
    *every* launch (`Chmy/src/KernelLaunch.jl:117`), letting consecutive kernels pipeline.

Safe because nothing ever reads a depth-integrated field off the `k = 1` plane: the 2D grid
operators (`∂x`, `∂y`, `lerp`/`hlerp` at `aa`/`ab`/`acx`/`acy`) never index `k ± 1`, and the
column kernels that read a 2D field index it explicitly — `_dotvel_staggered_bp!` writes
`lerp(H, NODE_ACX, grid2d, i, j, 1)`.
"""
struct FlatLauncher{WS,B,G}
    backend::B
    groupsize::G
    sync::Bool
end

function FlatLauncher(arch, grid; sync = true)
    bk = get_backend(arch)
    n = size(grid, Center())
    n[3] == 1 || throw(ArgumentError("FlatLauncher needs a depth-integrated grid, got nz = $(n[3])"))
    gs = heuristic_groupsize(bk, Val(3))
    return FlatLauncher{(n[1] + 2, n[2] + 2, 1),typeof(bk),typeof(gs)}(bk, gs, sync)
end

function (l::FlatLauncher{WS})(arch, grid, kernel_and_args::Pair; bc = nothing) where {WS}
    bc === nothing || throw(ArgumentError("FlatLauncher does not implement the bc= path"))
    kernel, args = kernel_and_args
    kernel(l.backend, l.groupsize, WS)(args..., Chmy.Offset(-1, -1, 0))
    l.sync && KernelAbstractions.synchronize(l.backend)
    return nothing
end

"""
    with_launcher(rt, launcher)

A copy of `rt` whose `launch2d` is `launcher`. `Runtime` holds its launchers as plain
fields, so a probe can swap the depth-integrated one without touching library code.
"""
with_launcher(rt, launcher) =
    Pagos.Runtime(rt.arch, rt.grid, rt.grid2d, rt.launch, launcher)

# ---------------------------------------------------------------------------
# Candidate kernels. Each is measured in `kernel_variants.jl` and stacked in `pt_loop.jl`;
# `pt_loop.jl` also checks that each reproduces `pseudo_transient!` bit-for-bit.
# ---------------------------------------------------------------------------

# Drops the two `copyto!`s in `update_basalstress!`: under the SSA limit `u_b = ū`, so this
# reads the depth-averaged velocity directly and writes `velocity.base_{x,y}` on the way
# past instead of copying into them first and reading them back.
@kernel inbounds = true function _basalstress_fused!(sbx, sby, bx, by, β, ux, uy, mask, grid, O)
    I = @index(Global, NTuple); I = I + O
    i, j, _ = I
    Z = zero(eltype(sbx))
    if node_active(mask, NODE_ACX, i, j)
        u = ux[I...]; bx[I...] = u; sbx[I...] = lerp(β, NODE_ACX, grid, I...) * u
    else
        bx[I...] = Z; sbx[I...] = Z
    end
    if node_active(mask, NODE_ACY, i, j)
        v = uy[I...]; by[I...] = v; sby[I...] = lerp(β, NODE_ACY, grid, I...) * v
    else
        by[I...] = Z; sby[I...] = Z
    end
end

# Drops the two `copyto!(u_old, u)`s: `u_old` is only needed as the pre-update iterate, so
# one kernel reads `u` once and writes both `u_old` and the relaxed `u`. Both components in
# one launch, halving the launch count as well.
@kernel inbounds = true function _vel_update_fused!(vx, vy, vx_old, vy_old,
                                                    dvx, dvy, dtx, dty, θ, O)
    I = @index(Global, NTuple); I = I + O
    a = vx[I...]; vx_old[I...] = a; vx[I...] = a + θ * dvx[I...] * dtx[I...]
    b = vy[I...]; vy_old[I...] = b; vy[I...] = b + θ * dvy[I...] * dty[I...]
end

# `_membrane_stress_staggered!` with the two viscosity/thickness prefactors read from
# precomputed fields. Both depend on `η` and `H` alone, which are fixed for the whole DIVA
# PT loop (`_iterate_viscosity!` is a no-op for DIVA; `NoDIVUpdate` leaves `µ̄` alone), so
# recomputing them per iteration recomputes `hlerp` — nine FP64 divisions a node — 200 times
# over for one answer.
@kernel inbounds = true function _membrane_pre!(sxx, sxy, syy, pre_aa, pre_ab, vel, mask, O)
    I = @index(Global, NTuple); I = I + O
    i, j, _ = I
    Z = zero(eltype(sxx))
    if node_active(mask, NODE_AA, i, j)
        p = pre_aa[I...]
        a = vel.depthaverage_x_dx[I...]; b = vel.depthaverage_y_dy[I...]
        sxx[I...] = p * (2a + b); syy[I...] = p * (a + 2b)
    else
        sxx[I...] = Z; syy[I...] = Z
    end
    sxy[I...] = node_fully_active(mask, NODE_AB, i, j) ?
        pre_ab[I...] * (vel.depthaverage_x_dy[I...] + vel.depthaverage_y_dx[I...]) : Z
end

@kernel inbounds = true function _membrane_prefactors!(pre_aa, pre_ab, η, H, mask, grid, O)
    I = @index(Global, NTuple); I = I + O
    i, j, _ = I
    Z = zero(eltype(pre_aa))
    pre_aa[I...] = node_active(mask, NODE_AA, i, j) ? 2 * η[I...] * H[I...] : Z
    # The kernel this replaces writes `2 * hlerp * lerp * (dxy + dyx) / 2`: the 2 and the
    # /2 cancel, so the prefactor is `hlerp * lerp`. Keeping the 2 here is a factor-of-two
    # error in `membrane_xy` that still converges — to a ~10%-different velocity field.
    pre_ab[I...] = node_fully_active(mask, NODE_AB, i, j) ?
        hlerp(η, NODE_AB, grid, I...) * lerp(H, NODE_AB, grid, I...) : Z
end

"""
    membrane_prefactor_fields(rt)

The `(pre_aa, pre_ab)` pair [`_membrane_pre!`](@ref) reads, at `aa` and `ab` on `grid2d`.
"""
membrane_prefactor_fields(rt) = (Field(rt.arch, rt.grid2d, Center()),
                                 Field(rt.arch, rt.grid2d, (Vertex(), Vertex(), Center())))

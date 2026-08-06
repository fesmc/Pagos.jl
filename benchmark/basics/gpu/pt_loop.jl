# ---------------------------------------------------------------------------
# Question: stacked on the real DIVA PT loop, what do the changes in
# `layout_and_launch.jl` and `kernel_variants.jl` add up to — and does the loop still
# produce the same answer?
#
# Four configurations, each timing a **fixed** 100 iterations (`abstol = 0`) so both sides
# do identical work and the comparison is throughput, not convergence:
#
#   baseline            `pseudo_transient!` as written, Chmy `Launcher`
#   + flat z-sweep      `FlatLauncher(sync = true)` swapped into `rt.launch2d`
#   + no per-launch sync `FlatLauncher(sync = false)`
#   + fused kernels     `pt_opt!` below — the same loop with three substitutions
#
# `pt_opt!` is a copy of `pseudo_transient!`'s `MomentumBalance2D` body with only those
# substitutions made; every other call (`drivingstress!`, `_tuning_init!`, `dotvel!`,
# `_tune!`, `bc!`, `_arm_tuning`, `_pt_error`) is the library's own, in the library's order.
# The last block checks each substitution reproduces `pseudo_transient!` **bit for bit**.
#
# !!! warning "Warm-up must reach the tuning cadence"
#     `AutotunedDynamicRelaxation(cadence = 50)` first calls `_arm_tuning`/`_tune!` at
#     iteration 50, so a 5-iteration warm-up leaves them — and their three-array
#     `mapreduce`s — uncompiled, and ~4 s of GPUCompiler then lands inside the *timed*
#     solve. That is the confound `docs/src/examples/ais-momentum/cpu-gpu.jl` hits.
#
# Run:  julia --project=benchmark benchmark/basics/gpu/pt_loop.jl [f32] [reps]
# ---------------------------------------------------------------------------
include(joinpath(@__DIR__, "common.jl"))

T    = length(ARGS) > 0 && ARGS[1] == "f32" ? Float32 : Float64
REPS = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 3

"""
    pt_opt!(mech, c, solver, rt, momentum, mask, pre_aa, pre_ab)

`pseudo_transient!` for a [`MomentumBalance2D`](@ref) with three substitutions:
[`_membrane_pre!`](@ref) (prefactors hoisted out of the loop), [`_basalstress_fused!`](@ref)
and [`_vel_update_fused!`](@ref). Valid for DIVA under `NoDIVUpdate`, where `η` and `H` do
not move — the `viscosity_continuation` hook is a no-op for DIVA by construction.
"""
function pt_opt!(mech, c, solver, rt, momentum, mask, pre_aa, pre_ab)
    (; velocity) = mech
    (; abstol, maxiter, ncheck, tuning) = solver
    ux, uy = velocity.depthaverage_x, velocity.depthaverage_y
    ux_old, uy_old = solver.velocity_x_old, solver.velocity_y_old
    st = mech.stress

    Pagos.drivingstress!(mech, c, rt, mask)
    state = Pagos._tuning_init!(tuning, solver, mech, c, rt, mask)
    scale = Pagos._convergence_scale(solver.convergence, mech, c, solver, rt, mask)
    rt.launch2d(rt.arch, rt.grid2d, _membrane_prefactors! =>
        (pre_aa, pre_ab, mech.material.viscosity_depthaveraged, mech.topography.thickness,
         mask, rt.grid2d))

    err = typemax(eltype(asarray(ux)))
    iter = 0
    while err > abstol && iter < maxiter
        iter += 1
        Pagos.depthaverage_velocitygradients!(velocity, rt, mask)
        rt.launch2d(rt.arch, rt.grid2d, _membrane_pre! =>
            (st.membrane_xx, st.membrane_xy, st.membrane_yy, pre_aa, pre_ab, velocity, mask))
        rt.launch2d(rt.arch, rt.grid2d, _basalstress_fused! =>
            (st.base_x, st.base_y, velocity.base_x, velocity.base_y,
             mech.friction.beta_eff, ux, uy, mask, rt.grid2d))
        Pagos.dotvel!(solver.velocity_x_dt, solver.velocity_y_dt, st.membrane_xx,
                      st.membrane_xy, st.membrane_yy, st.base_x, st.base_y, st.driving_x,
                      st.driving_y, mech.topography.thickness, c.density_ice, rt, momentum,
                      mask; gamma = state.gamma, resid_x = solver.residual_x,
                      resid_y = solver.residual_y)
        state = Pagos._tune!(tuning, state, solver, mech, c, rt, mask, ux, uy)
        rt.launch2d(rt.arch, rt.grid2d, _vel_update_fused! =>
            (ux, uy, ux_old, uy_old, solver.velocity_x_dt, solver.velocity_y_dt,
             solver.dtau_x, solver.dtau_y, state.theta_v))
        bc!(rt.arch, rt.grid2d, ux => Neumann())
        bc!(rt.arch, rt.grid2d, uy => Neumann())
        state = Pagos._arm_tuning(tuning, state, solver, ux, uy, iter)
        (iter % ncheck == 0 || iter == maxiter) &&
            (err = Pagos._pt_error(solver.convergence, solver, ux, uy, scale))
    end
    return (; iterations = iter, error = err, converged = err <= abstol)
end

# ---------------------------------------------------------------------------
const NITER = 100
fx  = gpu_fixture(; T)
mom = DIVAMomentumBalance()
pre_aa, pre_ab = membrane_prefactor_fields(fx.rt)

RT_BASE = fx.rt
RT_FLAT = with_launcher(fx.rt, FlatLauncher(fx.rt.arch, fx.rt.grid2d; sync = true))
RT_FAST = with_launcher(fx.rt, FlatLauncher(fx.rt.arch, fx.rt.grid2d; sync = false))

# One solver, reused. A fresh `PseudoTransientSolver` allocates eight `grid2d` fields (~22 MB
# each at this size), and building one per timed run exhausts the CUDA pool.
const SOLVER = PseudoTransientSolver(fx.grid; GPU_SOLVER_KWARGS..., maxiter = NITER,
                                     abstol = 0.0)
function reset_state!()
    for f in (SOLVER.velocity_x_dt, SOLVER.velocity_y_dt, SOLVER.residual_x,
              SOLVER.residual_y, SOLVER.velocity_x_old, SOLVER.velocity_y_old)
        setdata!(f, zero(T))
    end
    setdata!(fx.mech.velocity.depthaverage_x, zero(T))
    setdata!(fx.mech.velocity.depthaverage_y, zero(T))
    return nothing
end

function run_loop(rt, loop!)
    reset_state!()
    diva_update!(fx.mech, SOLVER, rt, fx.mask)
    CUDA.synchronize()
    t = CUDA.@elapsed loop!(rt)
    return 1e3 * t / NITER
end
base_loop(rt) = pseudo_transient!(fx.mech, fx.cst, SOLVER, rt, mom, fx.mask)
opt_loop(rt)  = pt_opt!(fx.mech, fx.cst, SOLVER, rt, mom, fx.mask, pre_aa, pre_ab)

VARIANTS = (("baseline (Chmy Launcher)",     () -> run_loop(RT_BASE, base_loop)),
            ("+ flat z-sweep",               () -> run_loop(RT_FLAT, base_loop)),
            ("+ no per-launch sync",         () -> run_loop(RT_FAST, base_loop)),
            ("+ fused kernels & prefactors", () -> run_loop(RT_FAST, opt_loop)))

for (nm, f) in VARIANTS            # compile everything off the clock
    print("compiling: ", nm, " ... "); flush(stdout)
    @printf("%.1f s\n", @elapsed f())
end
GC.gc(); CUDA.reclaim()

best = Dict{String,Float64}(); snap = Dict{String,Matrix{Float32}}()
for rep in 1:REPS, (label, f) in VARIANTS
    t = f()
    best[label] = min(get(best, label, Inf), t)
    rep == 1 && (snap[label] =
        Array(Float32.(asarray(fx.mech.velocity.depthaverage_x)[:, :, 1])))
end

@printf("\n=== %s, DIVA PT loop, %dx%d, %d fixed iterations (min of %d) ===\n",
        T, GPU_NX, GPU_NY, NITER, REPS)
b = best[first(VARIANTS)[1]]
for (label, _) in VARIANTS
    @printf("%-32s %7.3f ms/iter   %5.2fx\n", label, best[label], b / best[label])
end

ref = snap[first(VARIANTS)[1]]
println("\nvelocity field vs. baseline")
for (label, _) in VARIANTS
    @printf("  %-32s identical = %-6s  rel max|Δū| = %.3g\n", label,
            isequal(snap[label], ref),
            maximum(abs.(snap[label] .- ref)) / max(maximum(abs, ref), 1f-30))
end

#=

# DIVA: CPU vs. GPU

GPU-backed [`Chmy.Field`](@ref)s can't be filled element-by-element the way
[`run_solve`](@ref) fills CPU ones (`fill_from_grid!`/`fill_from_grid3d!`): CUDA.jl
disallows scalar indexing on a `CuArray`, since each write would become its own kernel
launch. Instead the whole [`MechanicState`](@ref) is built on the CPU as usual, then moved
to the GPU in one bulk `Adapt.adapt(CuArray, mech)` call — `Field` carries its own `Adapt`
rule for this, so the whole struct moves at once rather than field-by-field or
element-by-element.

[`IceMask`](@ref) wraps a `Field` reference rather than a copy, so it must be reconstructed
from the *adapted* topography's mask field rather than adapted itself.

A discarded warm-up solve runs first and is excluded from the timed comparison: GPU kernels
compile on first use (`GPUCompiler`), and that one-time cost is comparable in size to the
~20 s solves being measured. (`f32-f64.jl`'s Float64-vs-Float32 timing does not get this
treatment — that confound is flagged instead of fixed — because CPU JIT is a much smaller
fraction of a CPU solve's wall-clock.)

The warm-up runs **in place on `mech_gpu`** rather than on a second throwaway state: an 8 GB
card holds one adapted `MechanicState` (~3.3 GB) comfortably but not two at once, and
`GC.gc(); CUDA.reclaim()` on a discarded GPU state does not reliably return its memory to the
pool before a second allocation. Kernels are specialized on argument *type*, not identity, so
warming up on `mech_gpu` itself and resetting its velocity to zero before the timed solve
gets the same benefit without ever holding two full states at once.
=#
resolution_km = 8
include(joinpath(@__DIR__, "helpers.jl"))
using CUDA

diva = run_solve(DIVAMomentumBalance(), grid, rt, mask; SOLVER_KWARGS...)
println("DIVA, CPU: ", (; diva.converged, diva.iterations, diva.elapsed, diva.residual))

if CUDA.functional()
    arch_gpu = Arch(CUDABackend())
    grid_gpu = StaggeredGrid(arch_gpu, Float64, lx, ly, dx, dy, layering)
    rt_gpu   = Runtime(grid_gpu)

    mech_cpu_for_gpu = MechanicState(grid)
    fill_from_grid!(mech_cpu_for_gpu.topography.thickness, H_ice)
    fill_from_grid!(mech_cpu_for_gpu.topography.surface, z_srf)
    fill_from_grid!(mech_cpu_for_gpu.material.viscosity_depthaveraged, visc_bar)
    fill_from_grid!(mech_cpu_for_gpu.friction.beta_eff, beta_eff)
    fill_from_grid!(mech_cpu_for_gpu.friction.beta, beta)
    fill_from_grid3d!(mech_cpu_for_gpu.material.viscosity, visc3d)
    setdata!(mech_cpu_for_gpu.velocity.depthaverage_x, 0.0)
    setdata!(mech_cpu_for_gpu.velocity.depthaverage_y, 0.0)

    topo_gpu = Pagos.Adapt.adapt(CuArray, topo)
    mask_gpu = IceMask(topo_gpu.mask.is_momentum_solved)

    function solve_on!(mech, grid, rt, mask; solver_kwargs...)
        solver = PseudoTransientSolver(grid; solver_kwargs...)
        diva_update!(mech, solver, rt, mask)
        elapsed = @elapsed result = pseudo_transient!(mech, cst, solver, rt,
                                                       DIVAMomentumBalance(), mask)
        speed = Array(Float32.(sqrt.(
            (@views (interior(mech.velocity.depthaverage_x)[1:(end - 1), :, 1] .+
                     interior(mech.velocity.depthaverage_x)[2:end, :, 1]) ./ 2) .^ 2 .+
            (@views (interior(mech.velocity.depthaverage_y)[:, 1:(end - 1), 1] .+
                     interior(mech.velocity.depthaverage_y)[:, 2:end, 1]) ./ 2) .^ 2)))
        return (; speed, elapsed, iterations = result.iterations, converged = result.converged,
               residual = result.residual)
    end

    mech_gpu = Pagos.Adapt.adapt(CuArray, mech_cpu_for_gpu)
    mech_cpu_for_gpu = nothing
    GC.gc()
    CUDA.reclaim()

    # Warm-up: cheap throwaway solve (loose abstol, low maxiter) to force kernel compilation
    # off the clock — runs in place on mech_gpu (see note above); velocity resets to zero
    # before the real, timed solve reuses the same state.
    solve_on!(mech_gpu, grid_gpu, rt_gpu, mask_gpu; SOLVER_KWARGS..., maxiter = 5, abstol = 0.0)
    setdata!(mech_gpu.velocity.depthaverage_x, 0.0)
    setdata!(mech_gpu.velocity.depthaverage_y, 0.0)
    println("GPU warm-up done (kernels compiled, not timed).")

    diva_gpu = solve_on!(mech_gpu, grid_gpu, rt_gpu, mask_gpu; SOLVER_KWARGS...)
    mech_gpu = nothing
    GC.gc()
    println("DIVA, GPU: ", (; diva_gpu.converged, diva_gpu.iterations, diva_gpu.elapsed, diva_gpu.residual))

    diff_gpu = on_ice(diva.speed .- diva_gpu.speed)
    diff_gpu_vals = filter(!isnan, diff_gpu)
    @printf("CPU vs GPU DIVA (on-ice): RMSE = %.4g m/yr, max|Δ| = %.4g m/yr, elapsed %.3gs vs %.3gs (%.2gx)\n",
           sqrt(mean(abs2, diff_gpu_vals)), maximum(abs, diff_gpu_vals), diva.elapsed, diva_gpu.elapsed,
           diva.elapsed / diva_gpu.elapsed)

    fig4 = Figure(size = (1150, 620))
    crange4 = (0, quantile(filter(!isnan, on_ice(diva.speed)), 0.995))
    for (col, (data, title)) in enumerate(((on_ice(diva.speed), "DIVA, CPU"),
                                           (on_ice(diva_gpu.speed), "DIVA, GPU")))
        ax = Axis(fig4[1, col], xlabel = "x (km)", ylabel = col == 1 ? "y (km)" : "",
            aspect = DataAspect(), title = title)
        hm = heatmap!(ax, xc, yc, data; colorrange = crange4)
        col == 2 && Colorbar(fig4[1, 3], hm, label = "speed (m/yr)")
    end
    Label(fig4[2, 1:3],
        "CPU: $(diva.iterations) iter, $(round(diva.elapsed, digits = 2))s   |   " *
        "GPU: $(diva_gpu.iterations) iter, $(round(diva_gpu.elapsed, digits = 2))s   |   " *
        "speedup = $(round(diva.elapsed / diva_gpu.elapsed, sigdigits = 3))x   |   " *
        "RMSE = $(round(sqrt(mean(abs2, diff_gpu_vals)), sigdigits = 3)) m/yr",
        fontsize = 12)
    save("$figdir/cpu-gpu.png", fig4)
    fig4
else
    println("CUDA.functional() == false — skipping the CPU vs GPU comparison.")
end

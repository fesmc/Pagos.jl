# ---------------------------------------------------------------------------
# Question: with the launch worksize already fixed (`layout_and_launch.jl`), which kernel
# bodies on the DIVA/SSA PT residual path are still doing avoidable work?
#
# Four candidates, all measured at the flat worksize so the comparison is like-for-like:
#
#  0. Is `heuristic_groupsize`'s block shape actually the fastest one for the most
#     expensive kernel, or just a reasonable default nobody has measured against here?
#  1. `_membrane_stress_staggered!` recomputes `hlerp(η)` every iteration. Chmy's
#     `HarmonicLinear` rule is `inv(muladd(t, inv(b) - inv(a), inv(a)))` applied three times
#     for a 2D interpolation — **nine divisions a node**. On a consumer card FP64 runs at
#     1/32 rate, which is why this one kernel is the single most expensive thing in the
#     Float64 loop. Both the `aa` and the `ab` prefactor depend only on `η` and `H`, and
#     neither moves during a DIVA solve.
#  2. `update_basalstress!` copies `ū` into `velocity.base_{x,y}` and then reads it back.
#  3. `copyto!(u_old, u)` followed by `u = u_old + θ·dv·dτ` reads `u` twice.
#
# Run:  julia --project=benchmark benchmark/basics/gpu/kernel_variants.jl [f32]
# ---------------------------------------------------------------------------
include(joinpath(@__DIR__, "common.jl"))

T = length(ARGS) > 0 && ARGS[1] == "f32" ? Float32 : Float64
fx = gpu_fixture(; T)
(; mech, rt, mask) = fx
(; stress, velocity, material, topography, friction) = mech

bk = get_backend(rt.arch)
gs = heuristic_groupsize(bk, Val(3))
ws = (GPU_NX + 2, GPU_NY + 2, 1)
flat(k, args) = (k(bk, gs, ws)(args..., Chmy.Offset(-1, -1, 0)); nothing)

println("=== $T, $(GPU_NX)x$(GPU_NY) — all variants at the flat worksize $ws ===\n")

pre_aa, pre_ab = membrane_prefactor_fields(rt)
prefactor_args = (pre_aa, pre_ab, material.viscosity_depthaveraged, topography.thickness,
                  mask, rt.grid2d)
membrane_args = (stress.membrane_xx, stress.membrane_xy, stress.membrane_yy,
                 material.viscosity_depthaveraged, topography.thickness, velocity, mask,
                 rt.grid2d)
flat(_membrane_prefactors!, prefactor_args)

## ------------------------------------------------------------- 0. groupsize sweep
# `gs` above is `heuristic_groupsize`'s answer, never checked against alternatives on a
# real kernel. Swept here on the single most expensive kernel in the loop
# (`_membrane_stress_staggered!`, ~32% of device time, ~40% of roofline — see
# `README.md`'s "why" — so this is the kernel most worth tuning) at the same flat worksize
# every other section uses, so the comparison isolates groupsize alone.
#
# `g` is passed positionally, not read off a captured variable, so each candidate gets its
# own compiled `kernel(bk, g, ws)` — same reasoning as `src/api/runtime.jl`'s `FlatLauncher`
# taking `GroupSize` as a type parameter rather than a field (see its docstring): a
# `Val`/literal in the call signature reaches `StaticSize`, a value read off a box does not.
#
# !!! warning "This card needs round-robin timing, not `bench`'s consecutive batches"
#     `bench` times one candidate for `reps` batches, then moves to the next — on this
#     Max-Q RTX 2070 Super that measures clock drift, not groupsize. A confirmed thermal
#     soak (`nvidia-smi`: 930 → 1395 MHz, stable) did not fix it: the *same* `(32, 8, 1)`
#     config measured anywhere from 315 us to 595 us purely as a function of its position
#     in the candidate sequence — a 1.9x spread with the config held fixed. Every candidate
#     below instead gets one short timing per round, cycling through all candidates for
#     many rounds, so each accumulates measurements across the *same* drift trajectory and
#     drift cancels out of the comparison rather than aliasing onto whichever candidate
#     happens to run first.
println("(0) groupsize sweep on the flat worksize $ws — _membrane_stress_staggered!")
gs_candidates = ((32, 8, 1), (32, 4, 1), (32, 16, 1), (32, 32, 1), (16, 16, 1),
                 (16, 32, 1), (64, 4, 1), (64, 8, 1), (128, 2, 1), (128, 4, 1),
                 (256, 1, 1), (256, 2, 1))
flat_gs(k, args, g) = (k(bk, g, ws)(args..., Chmy.Offset(-1, -1, 0)); nothing)

for _ in 1:500; flat_gs(Pagos._membrane_stress_staggered!, membrane_args, gs); end
CUDA.synchronize()

n_per_round, nrounds = 10, 60
gs_samples = Dict(g => Float64[] for g in gs_candidates)
for _ in 1:nrounds, g in gs_candidates
    t = 1e6 * CUDA.@elapsed(begin
                                for _ in 1:n_per_round
                                    flat_gs(Pagos._membrane_stress_staggered!, membrane_args, g)
                                end
                            end) / n_per_round
    push!(gs_samples[g], t)
end
gs_medians = [g => median(gs_samples[g]) for g in gs_candidates]
for (g, med) in gs_medians
    @printf("  %-14s %6.1f us (median of %d rounds)\n", string(g), med, nrounds)
end
best_gs, best_t = gs_medians[argmin(last.(gs_medians))]
heur_t = gs_medians[findfirst(==(gs), first.(gs_medians))][2]
spread = maximum(last.(gs_medians)) / minimum(last.(gs_medians))
@printf("  -> heuristic %s: %.1f us; best %s: %.1f us; spread across all candidates %.2fx\n",
        gs, heur_t, best_gs, best_t, spread)
spread < 1.1 &&
    println("  -> within noise: this kernel is not groupsize-sensitive on this card.\n") ||
    println()

## ----------------------------------------------------- 1. membrane-stress prefactors
println("(1) membrane stress: hlerp per iteration vs. hoisted prefactors")
m1 = bench(() -> flat(Pagos._membrane_stress_staggered!, membrane_args),
           "_membrane_stress_staggered! (9 FP divisions)")
m2 = bench(() -> flat(_membrane_pre!, (stress.membrane_xx, stress.membrane_xy,
                                       stress.membrane_yy, pre_aa, pre_ab, velocity, mask)),
           "_membrane_pre! (prefactors read, 0 divisions)")
mp = bench(() -> flat(_membrane_prefactors!, prefactor_args),
           "  the prefactor build itself (once per solve)")
@printf("  -> %.2fx per iteration; the build pays for itself after %.1f iterations\n\n",
        m1 / m2, mp / (m1 - m2))

## ----------------------------------------------------------- 2. basal stress + copies
println("(2) update_basalstress!: two copyto!s + kernel vs. one fused kernel")
b1 = bench(() -> (copyto!(asarray(velocity.base_x), asarray(velocity.depthaverage_x));
                  copyto!(asarray(velocity.base_y), asarray(velocity.depthaverage_y));
                  flat(Pagos._basalstress_staggered!,
                       (stress.base_x, stress.base_y, friction.beta_eff, velocity.base_x,
                        velocity.base_y, mask, rt.grid2d))),
           "copyto! x2 + _basalstress_staggered!")
b2 = bench(() -> flat(_basalstress_fused!,
                      (stress.base_x, stress.base_y, velocity.base_x, velocity.base_y,
                       friction.beta_eff, velocity.depthaverage_x, velocity.depthaverage_y,
                       mask, rt.grid2d)), "_basalstress_fused!")
@printf("  -> %.2fx\n\n", b1 / b2)

## -------------------------------------------------------- 3. velocity update + copies
ux, uy = velocity.depthaverage_x, velocity.depthaverage_y
uxo = Field(rt.arch, rt.grid2d, NODE_ACX); uyo = Field(rt.arch, rt.grid2d, NODE_ACY)
dvx = Field(rt.arch, rt.grid2d, NODE_ACX); dvy = Field(rt.arch, rt.grid2d, NODE_ACY)
dtx = Field(rt.arch, rt.grid2d, NODE_ACX); dty = Field(rt.arch, rt.grid2d, NODE_ACY)
setdata!(dvx, T(1)); setdata!(dvy, T(1)); setdata!(dtx, T(1e-3)); setdata!(dty, T(1e-3))
θ = one(T)

println("(3) u_old <- u  then  u = u_old + θ·dv·dτ   (both components)")
s1 = bench(() -> begin
               copyto!(asarray(uxo), asarray(ux)); copyto!(asarray(uyo), asarray(uy))
               Pagos.pseudo_vel!(asarray(ux), asarray(uxo), asarray(dvx), asarray(dtx), θ)
               Pagos.pseudo_vel!(asarray(uy), asarray(uyo), asarray(dvy), asarray(dty), θ)
           end, "copyto! x2 + pseudo_vel! x2")
s2 = bench(() -> flat(_vel_update_fused!, (ux, uy, uxo, uyo, dvx, dvy, dtx, dty, θ)),
           "_vel_update_fused! (one launch)")
@printf("  -> %.2fx\n", s1 / s2)

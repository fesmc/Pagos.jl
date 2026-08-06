# ---------------------------------------------------------------------------
# Question: where does a depth-integrated Pagos kernel's GPU time actually go?
#
# Two answers that are properties of the *launch*, not of any kernel body:
#
#  1. A `Field` on `grid2d` is one k-plane of a `(nx + 4, ny + 4, 5)` parent — 5× the
#     memory its interior needs, because the halo convention adds two rings in z to an axis
#     with no extent. That is what puts a Float64 AIS state at 8 km over an 8 GB card.
#  2. Chmy's `Launcher` sweeps `size(grid, Center()) .+ 2`, i.e. `(nx + 2, ny + 2, 3)` on
#     `grid2d`. Every depth-integrated kernel therefore reads and writes **three** k-planes
#     where one would do.
#
# Run:  julia --project=benchmark benchmark/basics/gpu/layout_and_launch.jl [f32]
# ---------------------------------------------------------------------------
include(joinpath(@__DIR__, "common.jl"))

T = length(ARGS) > 0 && ARGS[1] == "f32" ? Float32 : Float64
fx = gpu_fixture(; T)
(; mech, rt, mask) = fx
(; stress, velocity, material, topography, friction) = mech

println("=== $T, $(GPU_NX)x$(GPU_NY), nz = $(GPU_NZ) ===\n")

## ------------------------------------------------------------------- field layout
println("field layout")
for (nm, f) in (("2D (viscosity_depthaveraged)", material.viscosity_depthaveraged),
                ("3D (viscosity)", material.viscosity))
    p, it = parent(f), interior(f)
    @printf("  %-30s interior %-16s parent %-16s %6.2f MB (%.1fx interior)\n",
            nm, string(size(it)), string(size(p)), sizeof(p) / 2^20, length(p) / length(it))
end
println("\n  launcher worksize  grid2d ", Chmy.KernelLaunch.worksize(rt.launch2d),
        "   grid ", Chmy.KernelLaunch.worksize(rt.launch))
println("  grid size (Center) grid2d ", size(rt.grid2d, Center()),
        "   grid ", size(rt.grid, Center()))

# Which k-planes a `grid2d` launch actually touches.
@kernel inbounds = true function _mark!(out, O)
    I = @index(Global, NTuple); I = I + O
    out[I...] = one(eltype(out))
end
scratch = Field(rt.arch, rt.grid2d, Center())
fill!(parent(scratch), zero(T))
rt.launch2d(rt.arch, rt.grid2d, _mark! => (scratch,))
pk = Array(parent(scratch))
println("  k-planes written by one grid2d launch: ",
        [k for k in axes(pk, 3) if any(!iszero, view(pk, :, :, k))], " of ", size(pk, 3))

## ---------------------------------------------------------------------- roofline
println("\ndevice bandwidth (the ceiling every kernel below is measured against)")
n2 = prod(size(parent(scratch))[1:2])
function bwline(label, f, bytes)
    us = bench(f; n = 100)
    @printf("  %-40s %8.1f us   %6.0f GB/s\n", label, us, bytes / (us * 1e-6) / 1e9)
end
a, b = CUDA.zeros(T, n2), CUDA.zeros(T, n2)
big, big2 = CUDA.zeros(T, 8n2), CUDA.zeros(T, 8n2)
bwline("contiguous copyto! (1 plane)", () -> copyto!(a, b), 2sizeof(a))
bwline("contiguous copyto! (8 planes)", () -> copyto!(big, big2), 2sizeof(big))
bwline("contiguous broadcast", () -> (a .= b .+ 1), 2sizeof(a))
big = big2 = nothing; GC.gc(); CUDA.reclaim()

ux, bx = velocity.depthaverage_x, velocity.base_x       # both at acx
iv_d, iv_s = interior(bx), interior(ux)
kp_d, kp_s = view(parent(bx), :, :, 3), view(parent(ux), :, :, 3)
bwline("copyto! on asarray(f) (strided view)", () -> copyto!(iv_d, iv_s),
       2sizeof(T) * length(iv_d))
bwline("copyto! on the parent k-plane", () -> copyto!(kp_d, kp_s), 2sizeof(T) * length(kp_d))

## ------------------------------------------------------- three k-planes versus one
bk = get_backend(rt.arch)
gs = heuristic_groupsize(bk, Val(3))
ws_full = size(rt.grid2d, Center()) .+ 2
ws_flat = (ws_full[1], ws_full[2], 1)
full(k, args) = (k(bk, gs, ws_full)(args..., Chmy.Offset(-1)); nothing)
flat(k, args) = (k(bk, gs, ws_flat)(args..., Chmy.Offset(-1, -1, 0)); nothing)

KERNELS = (
    ("_depthaverage_velocity_gradients!", Pagos._depthaverage_velocity_gradients!,
     (velocity, mask, rt.grid2d)),
    ("_membrane_stress_staggered!", Pagos._membrane_stress_staggered!,
     (stress.membrane_xx, stress.membrane_xy, stress.membrane_yy,
      material.viscosity_depthaveraged, topography.thickness, velocity, mask, rt.grid2d)),
    ("_basalstress_staggered!", Pagos._basalstress_staggered!,
     (stress.base_x, stress.base_y, friction.beta_eff, velocity.base_x, velocity.base_y,
      mask, rt.grid2d)),
)

println("\nthe same kernel over 3 k-planes (Chmy default) and over 1")
tot_full = tot_flat = 0.0
for (nm, k, args) in KERNELS
    tf = bench(() -> full(k, args), "$nm  [3 planes]")
    tl = bench(() -> flat(k, args), "$nm  [1 plane]")
    @printf("  %-46s %8.2fx\n", "->", tf / tl)
    global tot_full += tf; global tot_flat += tl
end
@printf("\n  hot depth-integrated kernels: %.0f us -> %.0f us  (%.2fx)\n",
        tot_full, tot_flat, tot_full / tot_flat)

# The flat launch must leave the interior untouched — the halo planes are the only
# difference, and nothing reads them.
Pagos.membranestress!(stress, velocity, material, topography, DIVAMomentumBalance(), rt, mask)
ref = copy(Array(asarray(stress.membrane_xy)))
fill!(parent(stress.membrane_xy), T(NaN))
flat(Pagos._membrane_stress_staggered!, KERNELS[2][3])
CUDA.synchronize()
println("\n  flat launch reproduces the full launch's interior: ",
        isequal(ref, Array(asarray(stress.membrane_xy))))

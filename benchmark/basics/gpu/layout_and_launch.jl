# ---------------------------------------------------------------------------
# Question: where does a depth-integrated Pagos kernel's GPU time actually go?
#
# Two answers that were properties of the *launch*, not of any kernel body:
#
#  1. A `Field` on `grid2d` used to be one k-plane of a `(nx + 4, ny + 4, 5)` parent — 5× the
#     memory its interior needs, because the halo convention adds two rings in z to an axis
#     with no extent. That is what put a Float64 AIS state at 8 km over an 8 GB card.
#  2. Chmy's `Launcher` sweeps `size(grid, Center()) .+ 2`, i.e. `(nx + 2, ny + 2, 3)` on
#     `grid2d`. Every depth-integrated kernel therefore read and wrote **three** k-planes
#     where one would do.
#
# **Both landed in `src`** (`pagos-roadmap/memreduce.md`): `rt.launch2d` is a
# `Pagos.FlatLauncher` (finding 2), and `MechanicState`'s `grid2d` fields now allocate with
# `halo = (h, h, 0)` (finding 1, via `Pagos._halo2d`). Measured numbers for both live in
# `benchmark/basics/gpu/README.md` §1, taken while getting the fix right — not reproduced
# live below, because that would mean reconstructing the pre-fix 5-plane allocation just to
# re-demonstrate a bug `mech`'s own fields can no longer exhibit (a `full()` launch below
# would now be a `BoundsError`, correctly: `mech`'s `grid2d` fields have no `z ± 1` to read).
# What is still live: the plane-count probe (a standalone `halo = 1` field, so it keeps
# showing what an *ordinary* Chmy `Launcher` sweeps) and the roofline comparisons, neither of
# which depended on the bug.
#
# Run:  julia --project=benchmark benchmark/basics/gpu/layout_and_launch.jl [f32]
# ---------------------------------------------------------------------------
include(joinpath(@__DIR__, "common.jl"))

T = length(ARGS) > 0 && ARGS[1] == "f32" ? Float32 : Float64
fx = gpu_fixture(; T)
(; mech, rt, mask) = fx
(; stress, velocity, material, friction) = mech

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
# Still a 5-plane parent, unlike most of `mech`'s other `grid2d` fields: `depthaverage_x`
# and `base_x` are 2 of the 14 fields `MechanicState` exempts from the z-ghost shrink
# because they're `bc!`'d (see its constructor's note) — so the data plane is still at
# parent index 3, exactly as a `halo = 1` field puts it.
kp_d, kp_s = view(parent(bx), :, :, 3), view(parent(ux), :, :, 3)
bwline("copyto! on asarray(f) (strided view)", () -> copyto!(iv_d, iv_s),
       2sizeof(T) * length(iv_d))
bwline("copyto! on the parent k-plane", () -> copyto!(kp_d, kp_s), 2sizeof(T) * length(kp_d))

## ------------------------------------------------ finding 1 is no longer reproducible
# What used to be here timed the same kernel launched over 3 k-planes (Chmy's default
# `Launcher` worksize) against 1 (`FlatLauncher`'s). That comparison needed a `grid2d` field
# with a real z-ghost to sweep the extra two planes of — `mech`'s fields no longer have one,
# so the `full()` launch below is now a `BoundsError`, on purpose: it is the fix, observed
# directly rather than timed. The historical numbers (2.9–4.2× per kernel) are in
# `benchmark/basics/gpu/README.md` §1.
bk = get_backend(rt.arch)
gs = heuristic_groupsize(bk, Val(3))
ws_full = size(rt.grid2d, Center()) .+ 2
full(k, args) = (k(bk, gs, ws_full)(args..., Chmy.Offset(-1)); nothing)
try
    full(Pagos._basalstress_staggered!,
         (stress.base_x, stress.base_y, friction.beta_eff, velocity.base_x, velocity.base_y,
          mask, rt.grid2d))
    println("\nunexpected: a 3-plane launch over a grid2d field succeeded")
catch e
    println("\na 3-plane launch over a grid2d field now errors, as expected: ", nameof(typeof(e)))
end

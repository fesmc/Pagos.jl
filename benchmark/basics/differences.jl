using Pagos
using Chairmarks
using Printf

# ---------------------------------------------------------------------------
# Benchmark: ∂x₁₂ (fused) vs sequential ∂x₁ + ∂x₂
#
# These kernels are memory-bandwidth-bound (~1 FLOP/byte).  The fused kernel
# reads `u` only once rather than twice, so it should be roughly 2× faster
# when `u` does not fit in cache.
#
# CPU-only benchmark: KA CPU kernels complete before returning, so no explicit
# synchronization is needed in the timed expressions. For GPU timings see
# difference_kernels.jl, which synchronizes inside the timed region (operators
# launch asynchronously on GPU).
# ---------------------------------------------------------------------------

const SIZES = [(128, 128), (256, 256), (512, 512), (1024, 1024)]

w = 72
println("=" ^ w)
println(" ∂x₁ / ∂x₂ / ∂x₁₂ benchmark — CPU (KernelAbstractions CPU backend)")
println("=" ^ w)
@printf("  %-14s %11s %11s %11s %9s\n",
    "Grid", "∂x₁ [μs]", "∂x₂ [μs]", "∂x₁₂ [μs]", "seq/fused")
println("-" ^ w)

for (nx, ny) in SIZES
    u   = rand(Float64, nx, ny)
    du₁ = similar(u)
    du₂ = similar(u)

    idx₁ = FlatIndexing(1, nx)
    idx₂ = FlatIndexing(1, ny)
    t1  = (@b ∂x!($du₁, $u, 1.0, $idx₁)).time
    t2  = (@b ∂y!($du₂, $u, 1.0, $idx₂)).time
    t12 = (@b ∂x₁₂!($du₁, $du₂, $u, 1.0, 1.0, $idx₁, $idx₂)).time

    @printf("  %-14s %11.1f %11.1f %11.1f %9.2fx\n",
        "$(nx)×$(ny)",
        t1  * 1e6,
        t2  * 1e6,
        t12 * 1e6,
        (t1 + t2) / t12)
end

println("=" ^ w)
println()

# ---------------------------------------------------------------------------
# Benchmark: indexing types — overhead of different boundary conditions
# ---------------------------------------------------------------------------

println("=" ^ w)
println(" ∂x₁ benchmark — indexing types ($(SIZES[end][1])×$(SIZES[end][2]) grid)")
println("=" ^ w)
@printf "  %-28s %11s\n" "Indexing" "∂x₁ [μs]"
println("-" ^ w)

let (nx, ny) = SIZES[end]
    u   = rand(Float64, nx, ny)
    du  = similar(u)
    idxs = [
        ("FlatIndexing",        FlatIndexing(1, nx)),
        ("ReflectiveIndexing",  ReflectiveIndexing(1, nx)),
        ("PeriodicIndexing",    PeriodicIndexing(1, nx)),
        ("StrictIndexing",      StrictIndexing(1, nx)),
    ]
    for (name, idx) in idxs
        t = (@b ∂x!($du, $u, 1.0, $idx)).time
        @printf "  %-28s %11.1f\n" name t * 1e6
    end
end

println("=" ^ w)
println()
println("Notes:")
println("  ∂x₁₂  reads u once; ∂x₁ + ∂x₂ reads u twice — fused wins when u > L3 cache.")
println("  All indexing variants compile to a single branch-free stencil via stencil_fd.")

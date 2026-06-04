using Pagos
using Chairmarks
using Printf
using CUDA

const HAS_CUDA = CUDA.functional()

# ---------------------------------------------------------------------------
# isbits functors used across all benchmarks
# ---------------------------------------------------------------------------

struct Addition end
g(x, y, ::Addition) = x + y    # dispatch-param style
add = Addition()

struct GAdd end
(::GAdd)(x, y) = x + y         # functor style (no scalar param)
gadd = GAdd()

# ---------------------------------------------------------------------------
# Benchmark 1: apply!(functor) vs full-array broadcasting
#
# For a mask that covers fraction f of the grid, apply! does ~f×N work while
# broadcasting does N work.  The cross-over depends on the mask density and
# kernel complexity.  Here the mask covers the central 50% square (~25% of
# all cells).
# ---------------------------------------------------------------------------

const SIZES = [(256, 256), (512, 512), (1024, 1024)]

w = 68
println("=" ^ w)
println(" apply! vs broadcasting — CPU, functor-dispatch (25 % active cells)")
println("=" ^ w)
@printf "  %-14s %14s %14s %10s\n" "Grid" "apply! [μs]" "broadcast [μs]" "ratio"
println("-" ^ w)

for (nx, ny) in SIZES
    x   = rand(Float32, nx, ny)
    y   = rand(Float32, nx, ny)
    out = zeros(Float32, nx, ny)

    mask = falses(nx, ny)
    mask[nx÷4:3nx÷4, ny÷4:3ny÷4] .= true
    acm  = ActiveCellsMap(mask)

    t_apply = (@b apply!($gadd, $acm, $out, ($x, $y))).time
    t_bcast = (@b $out .= $x .+ $y).time

    @printf("  %-14s %14.1f %14.1f %9.2fx\n",
        "$(nx)×$(ny)", t_apply * 1e6, t_bcast * 1e6, t_bcast / t_apply)
end

println("=" ^ w)
println()

# ---------------------------------------------------------------------------
# Benchmark 2: apply! dispatch styles — functor vs dispatch-param
# ---------------------------------------------------------------------------

println("=" ^ w)
println(" apply! dispatch styles — CPU ($(SIZES[end][1])×$(SIZES[end][2]), 25 % active)")
println("=" ^ w)
@printf "  %-32s %14s\n" "Style" "time [μs]"
println("-" ^ w)

let (nx, ny) = SIZES[end]
    x   = rand(Float32, nx, ny)
    y   = rand(Float32, nx, ny)
    out = zeros(Float32, nx, ny)
    mask = falses(nx, ny)
    mask[nx÷4:3nx÷4, ny÷4:3ny÷4] .= true
    acm  = ActiveCellsMap(mask)

    t1 = (@b apply!($gadd, $acm, $out, ($x, $y))).time
    t2 = (@b apply!($g,    $acm, $out, ($x, $y), ($add,))).time
    t3 = (@b $out .= $x .+ $y).time

    @printf "  %-32s %14.1f\n" "apply!(GAdd())" t1 * 1e6
    @printf "  %-32s %14.1f\n" "apply!(g, ..., (Addition(),))" t2 * 1e6
    @printf "  %-32s %14.1f\n" "out .= x .+ y  (full array)" t3 * 1e6
end

println("=" ^ w)
println()

# ---------------------------------------------------------------------------
# Benchmark 3: apply! vs creep! vs broadcast for SmithMorlandCreep
# ---------------------------------------------------------------------------

println("=" ^ w)
println(" SmithMorlandCreep — apply! vs creep! vs broadcast — CPU")
println("=" ^ w)
@printf "  %-14s %14s %14s %14s\n" "Grid" "apply! [μs]" "creep! [μs]" "bcast [μs]"
println("-" ^ w)

smc = SmithMorlandCreep()
(; p0, p2, p4, D_0, σ_0) = smc

for (nx, ny) in SIZES
    x   = rand(Float32, nx, ny)
    out = zeros(Float32, nx, ny)
    mask = falses(nx, ny)
    mask[nx÷4:3nx÷4, ny÷4:3ny÷4] .= true
    acm  = ActiveCellsMap(mask)

    t_apply = (@b apply!(creep, $acm, $out, ($x,), ($smc,))).time
    t_creep = (@b creep!($out, $x, $smc)).time
    t_bcast = (@b $out .= $D_0 ./ $x .* ($p0 .+ $p2 .* ($x ./ $σ_0).^2 .+ $p4 .* ($x ./ $σ_0).^4)).time

    @printf("  %-14s %14.1f %14.1f %14.1f\n",
        "$(nx)×$(ny)", t_apply * 1e6, t_creep * 1e6, t_bcast * 1e6)
end

println("=" ^ w)
println()

# ---------------------------------------------------------------------------
# GPU benchmarks (skipped if no CUDA device is available)
# ---------------------------------------------------------------------------

if HAS_CUDA
    println("=" ^ w)
    println(" apply! vs broadcasting — GPU (25 % active cells)")
    println("=" ^ w)
    @printf "  %-14s %14s %14s %10s\n" "Grid" "apply! [μs]" "broadcast [μs]" "ratio"
    println("-" ^ w)

    for (nx, ny) in SIZES
        x   = CUDA.rand(Float32, nx, ny)
        y   = CUDA.rand(Float32, nx, ny)
        out = CUDA.zeros(Float32, nx, ny)

        mask_cpu = falses(nx, ny)
        mask_cpu[nx÷4:3nx÷4, ny÷4:3ny÷4] .= true
        acm = ActiveCellsMap(cu(mask_cpu))

        t_apply = (@b begin
            apply!($gadd, $acm, $out, ($x, $y))
            CUDA.synchronize()
        end).time
        t_bcast = (@b begin
            $out .= $x .+ $y
            CUDA.synchronize()
        end).time

        @printf("  %-14s %14.1f %14.1f %9.2fx\n",
            "$(nx)×$(ny)", t_apply * 1e6, t_bcast * 1e6, t_bcast / t_apply)
    end

    println("=" ^ w)
    println()

    println("=" ^ w)
    println(" SmithMorlandCreep — apply! vs creep! vs broadcast — GPU")
    println("=" ^ w)
    @printf "  %-14s %14s %14s %14s\n" "Grid" "apply! [μs]" "creep! [μs]" "bcast [μs]"
    println("-" ^ w)

    for (nx, ny) in SIZES
        x   = CUDA.rand(Float32, nx, ny)
        out = CUDA.zeros(Float32, nx, ny)
        mask_cpu = falses(nx, ny)
        mask_cpu[nx÷4:3nx÷4, ny÷4:3ny÷4] .= true
        acm = ActiveCellsMap(cu(mask_cpu))

        t_apply = (@b begin
            apply!(creep, $acm, $out, ($x,), ($smc,))
            CUDA.synchronize()
        end).time
        t_creep = (@b begin
            creep!($out, $x, $smc)
            CUDA.synchronize()
        end).time
        t_bcast = (@b begin
            $out .= $D_0 ./ $x .* ($p0 .+ $p2 .* ($x ./ $σ_0).^2 .+ $p4 .* ($x ./ $σ_0).^4)
            CUDA.synchronize()
        end).time

        @printf("  %-14s %14.1f %14.1f %14.1f\n",
            "$(nx)×$(ny)", t_apply * 1e6, t_creep * 1e6, t_bcast * 1e6)
    end

    println("=" ^ w)
else
    println("GPU benchmarks skipped: no CUDA-capable device found (CUDA.functional() == false)")
end

println()
println("Notes:")
println("  apply!  operates only on active (masked) cells — ~25 % of the grid here.")
println("  creep!  operates on the full array; apply! wins when the mask is sparse.")
println("  bcast   is the unmasked baseline: reads and writes every element.")

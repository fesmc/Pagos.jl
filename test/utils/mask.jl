using Test
using Pagos
using CUDA

# Init struct for isbit functor
struct ScaleKernel{T}
    factor::T
end
(k::ScaleKernel)(x) = k.factor * x

# Init struct for dispatch test
struct Addition end
struct Subtraction end

# Functions for dispatch test
g(x, y, ::Addition)    = x + y
g(x, y, ::Subtraction) = x - y

@testset "ActiveCellsMap construction" begin
    mask = Bool[1 0 1; 0 1 0; 1 1 0]
    acm  = ActiveCellsMap(mask)

    @test length(acm.indices) == count(mask)
    @test Set(acm.indices) == Set(findall(mask))
end

@testset "apply! single input" begin
    mask = Bool[1 0 1; 0 1 0; 1 1 0]
    x    = Float64[1 2 3; 4 5 6; 7 8 9]
    out  = zeros(Float64, size(x))
    acm  = ActiveCellsMap(mask)

    apply!(ScaleKernel(2.0), acm, out, x)

    for I in CartesianIndices(mask)
        if mask[I]
            @test out[I] ≈ 2.0 * x[I]
        else
            @test out[I] == 0.0
        end
    end
end

@testset "apply! two inputs" begin
    mask = Bool[1 0 1; 0 1 0; 1 1 0]
    x    = Float64[1 2 3; 4 5 6; 7 8 9]
    y    = Float64[9 8 7; 6 5 4; 3 2 1]
    out  = zeros(Float64, size(x))
    acm  = ActiveCellsMap(mask)

    apply!(+, acm, out, x, y)

    for I in CartesianIndices(mask)
        if mask[I]
            @test out[I] ≈ x[I] + y[I]
        else
            @test out[I] == 0.0
        end
    end
end

@testset "apply! empty mask" begin
    mask = falses(4, 4)
    out  = ones(Float64, 4, 4)
    acm  = ActiveCellsMap(mask)

    apply!(ScaleKernel(99.0), acm, out, out)

    @test all(out .== 1.0)
end

@testset "apply! with dispatch param" begin
    mask = Bool[1 0 1; 0 1 0; 1 1 0]
    x    = Float64[1 2 3; 4 5 6; 7 8 9]
    y    = Float64[9 8 7; 6 5 4; 3 2 1]
    C    = zeros(Float64, size(x))
    D    = zeros(Float64, size(x))
    acm  = ActiveCellsMap(mask)

    apply!(g, acm, C, (x, y), (Addition(),))
    apply!(g, acm, D, (x, y), (Subtraction(),))

    for I in CartesianIndices(mask)
        if mask[I]
            @test C[I] ≈ x[I] + y[I]
            @test D[I] ≈ x[I] - y[I]
        else
            @test C[I] == 0.0
            @test D[I] == 0.0
        end
    end
end

@testset "apply! with dispatch and cuarrays" begin
    mask = cu(Bool[1 0 1; 0 1 0; 1 1 0])
    x    = cu(Float32[1 2 3; 4 5 6; 7 8 9])
    y    = cu(Float32[9 8 7; 6 5 4; 3 2 1])
    C    = CUDA.zeros(Float32, size(x))
    D    = CUDA.zeros(Float32, size(x))
    acm  = ActiveCellsMap(mask)

    apply!(g, acm, C, (x, y), (Addition(),))
    apply!(g, acm, D, (x, y), (Subtraction(),))

    @test count(C .== x .+ y) == count(mask)
    @test count(D .== x .- y) == count(mask)
end

x = CUDA.rand(1000, 1000)
y = CUDA.rand(1000, 1000)
C = CUDA.zeros(Float32, size(x))
mask = CUDA.fill(false, size(x))
mask[250:750, 250:750] .= true
acm = ActiveCellsMap(mask)

g(x, y, p::Addition) = x + y
add = Addition()

struct GAdd end
(gadd::GAdd)(x, y) = x + y
gadd = GAdd()

@b apply!($g, $acm, $C, ($x, $y), ($add,))
@b apply!($gadd, $acm, $C, ($x, $y))
@b $x .+ $y

# This really show that map provides the most efficient implementation, both in terms of computation time and memeory allocation
smc = SmithMorlandCreep()
(; p0, p2, p4, D_0, σ_0) = smc
@b apply!(creep, $acm, $C, ($x,), (smc,))
@b $D_0 ./ $x .* ($p0 .+ $p2 .* ($x ./ $σ_0).^2 .+ $p4 .* ($x ./ $σ_0).^4)
@b creep!($C, $x, smc)

field_1km = rand(Float32, 6000, 6000)
Base.format_bytes(sizeof(field_1km))
100 * 137   # 100 fields ==> about 14 GB, still fits on GPU memory
using Pagos
using Test

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Linear field u[i,j] = float(i), exact ∂u/∂x₁ = 1/dx everywhere.
linear_x1(nx, ny) = Float64[float(i) for i in 1:nx, _ in 1:ny]
linear_x2(nx, ny) = Float64[float(j) for _ in 1:nx, j in 1:ny]
linear_x3(nx, ny, nz) = Float64[float(k) for _ in 1:nx, _ in 1:ny, k in 1:nz]

# ---------------------------------------------------------------------------
# ∂x₁ — first-dimension derivative
# ---------------------------------------------------------------------------

@testset "∂x₁ constant field → zero" begin
    u  = ones(Float64, 8, 6)
    du = ∂x₁(u, 0.5)
    @test all(du .≈ 0.0)
end

@testset "∂x₁ linear field — FlatIndexing (exact everywhere)" begin
    # FlatIndexing clamps: one-sided h=1 at boundaries, central h=2 interior.
    # For u[i,j]=i both give du = 1/dx exactly.
    nx, ny = 8, 6
    dx = 0.5
    u  = linear_x1(nx, ny)
    du = ∂x₁(u, dx)
    @test all(du .≈ 1.0 / dx)
end

@testset "∂x₁ of x₂-only field → zero" begin
    nx, ny = 8, 6
    u  = linear_x2(nx, ny)
    du = ∂x₁(u, 1.0)
    @test all(du .≈ 0.0)
end

@testset "∂x₁ ReflectiveIndexing — zero gradient at boundary" begin
    # ReflectiveIndexing mirrors the grid: at i=1, im1 = ip1 = 2 → (u[2]-u[2])/(2dx) = 0.
    nx, ny = 8, 6
    dx = 1.0
    u  = linear_x1(nx, ny)
    du = ∂x₁(u, dx, ReflectiveIndexing(1, nx))
    @test all(du[1,   :] .≈ 0.0)
    @test all(du[end, :] .≈ 0.0)
    # Interior is still exact central difference.
    @test all(du[2:end-1, :] .≈ 1.0 / dx)
end

@testset "∂x₁ PeriodicIndexing — interior exact, boundaries wrap" begin
    nx, ny = 8, 6
    dx = 1.0
    u  = linear_x1(nx, ny)
    du = ∂x₁(u, dx, PeriodicIndexing(1, nx))
    # Interior cells see a symmetric two-cell stencil from the same linear field.
    @test all(du[2:end-1, :] .≈ 1.0 / dx)
    # Boundary cells wrap: i=1 → im1=nx, i=nx → ip1=1; not 1/dx for a non-periodic field.
    @test !all(du[1, :] .≈ 1.0 / dx)
end

@testset "∂x₁ in-place matches allocating" begin
    nx, ny = 10, 7
    u   = rand(Float64, nx, ny)
    du1 = ∂x₁(u, 0.4)
    du2 = similar(u)
    ∂x₁!(du2, u, 0.4)
    @test du1 ≈ du2
end

# ---------------------------------------------------------------------------
# ∂x₂ — second-dimension derivative
# ---------------------------------------------------------------------------

@testset "∂x₂ constant field → zero" begin
    u  = ones(Float64, 8, 6)
    du = ∂x₂(u, 0.3)
    @test all(du .≈ 0.0)
end

@testset "∂x₂ linear field — FlatIndexing (exact everywhere)" begin
    nx, ny = 8, 6
    dy = 0.25
    u  = linear_x2(nx, ny)
    du = ∂x₂(u, dy)
    @test all(du .≈ 1.0 / dy)
end

@testset "∂x₂ of x₁-only field → zero" begin
    nx, ny = 8, 6
    u  = linear_x1(nx, ny)
    du = ∂x₂(u, 1.0)
    @test all(du .≈ 0.0)
end

@testset "∂x₂ ReflectiveIndexing — zero gradient at boundary" begin
    nx, ny = 8, 6
    dy = 1.0
    u  = linear_x2(nx, ny)
    du = ∂x₂(u, dy, ReflectiveIndexing(1, ny))
    @test all(du[:,   1] .≈ 0.0)
    @test all(du[:, end] .≈ 0.0)
    @test all(du[:, 2:end-1] .≈ 1.0 / dy)
end

@testset "∂x₂ in-place matches allocating" begin
    nx, ny = 10, 7
    u   = rand(Float64, nx, ny)
    du1 = ∂x₂(u, 0.6)
    du2 = similar(u)
    ∂x₂!(du2, u, 0.6)
    @test du1 ≈ du2
end

# ---------------------------------------------------------------------------
# ∂x₃ — third-dimension derivative
# ---------------------------------------------------------------------------

@testset "∂x₃ linear field — FlatIndexing (exact everywhere)" begin
    nx, ny, nz = 4, 5, 6
    dz = 2.0
    u  = linear_x3(nx, ny, nz)
    du = ∂x₃(u, dz)
    @test all(du .≈ 1.0 / dz)
end

@testset "∂x₃ in-place matches allocating" begin
    nx, ny, nz = 4, 5, 6
    u   = rand(Float64, nx, ny, nz)
    du1 = ∂x₃(u, 1.5)
    du2 = similar(u)
    ∂x₃!(du2, u, 1.5)
    @test du1 ≈ du2
end

# ---------------------------------------------------------------------------
# ∂x₃ with sigma transform
# ---------------------------------------------------------------------------

# For any field linear in sigma (u[i,j,k] = ζ_aa[k]), central differences in sigma
# space are exact and the physical derivative is ∂u/∂z = (∂u/∂ζ)/H = 1/H.
# This holds at every grid point (boundary and interior) regardless of the
# sigma distribution, so it is the canonical test for the sigma dispatch.

@testset "∂x₃ sigma — constant field → zero (PowerSigmaTransform)" begin
    nx, ny, nz = 4, 5, 8
    transform = PowerSigmaTransform{Float64}(nz, 2)
    H  = rand(Float64, nx, ny) .+ 100.0
    u  = ones(Float64, nx, ny, nz)
    du = ∂x₃(u, H, transform)
    @test all(du .≈ 0.0)
end

@testset "∂x₃ sigma — linear-in-sigma → exact 1/H (PowerSigmaTransform)" begin
    nx, ny, nz = 4, 5, 8
    transform = PowerSigmaTransform{Float64}(nz, 2)
    ζ = get_ζ_aa(Float64, transform)
    H = rand(Float64, nx, ny) .+ 100.0
    u = Float64[ζ[k] for _ in 1:nx, _ in 1:ny, k in 1:nz]
    du = ∂x₃(u, H, transform)
    expected = Float64[1.0 / H[i, j] for i in 1:nx, j in 1:ny, _ in 1:nz]
    @test du ≈ expected
end

@testset "∂x₃ sigma in-place matches allocating" begin
    nx, ny, nz = 4, 5, 8
    transform = PowerSigmaTransform{Float64}(nz, 2)
    H   = rand(Float64, nx, ny) .+ 100.0
    u   = rand(Float64, nx, ny, nz)
    du1 = ∂x₃(u, H, transform)
    du2 = similar(u)
    ∂x₃!(du2, u, H, transform)
    @test du1 ≈ du2
end

@testset "∂x₃ sigma agrees with uniform ∂x₃ when sigma is linear and H is constant" begin
    # PowerSigmaTransform with exponent=1 gives uniform sigma spacing dζ = 1/(nz-1).
    # With constant thickness H₀, ∂u/∂z via sigma == ∂u/∂z via uniform dz = H₀/(nz-1).
    nx, ny, nz = 4, 5, 8
    H₀ = 500.0
    dz = H₀ / (nz - 1)
    transform = PowerSigmaTransform{Float64}(nz, 1)   # linear (uniform) sigma
    H  = fill(H₀, nx, ny)
    u  = rand(Float64, nx, ny, nz)
    du_sigma   = ∂x₃(u, H, transform)
    du_uniform = ∂x₃(u, dz)
    @test du_sigma ≈ du_uniform
end

# ---------------------------------------------------------------------------
# ∂x₁₂ — fused planar derivative
# ---------------------------------------------------------------------------

@testset "∂x₁₂ matches separate ∂x₁ and ∂x₂" begin
    nx, ny = 10, 8
    dx, dy = 0.4, 0.6
    u = rand(Float64, nx, ny)
    du₁_sep = ∂x₁(u, dx)
    du₂_sep = ∂x₂(u, dy)
    du₁_fus, du₂_fus = ∂x₁₂(u, dx, dy)
    @test du₁_fus ≈ du₁_sep
    @test du₂_fus ≈ du₂_sep
end

@testset "∂x₁₂! in-place matches allocating" begin
    nx, ny = 10, 8
    dx, dy = 0.4, 0.6
    u   = rand(Float64, nx, ny)
    du₁_a, du₂_a = ∂x₁₂(u, dx, dy)
    du₁_b = similar(u)
    du₂_b = similar(u)
    ∂x₁₂!(du₁_b, du₂_b, u, dx, dy)
    @test du₁_b ≈ du₁_a
    @test du₂_b ≈ du₂_a
end

@testset "∂x₁₂ per-dimension indexing" begin
    nx, ny = 12, 10
    dx, dy = 1.0, 1.0
    u = rand(Float64, nx, ny)
    # Mix indexing types for each dimension independently.
    idx₁ = PeriodicIndexing(1, nx)
    idx₂ = ReflectiveIndexing(1, ny)
    du₁_fus, du₂_fus = ∂x₁₂(u, dx, dy, idx₁, idx₂)
    du₁_sep = ∂x₁(u, dx, idx₁)
    du₂_sep = ∂x₂(u, dy, idx₂)
    @test du₁_fus ≈ du₁_sep
    @test du₂_fus ≈ du₂_sep
end
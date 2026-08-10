using Pagos
using Test
using KernelAbstractions: @kernel, @index

###############################################################
# Manufactured solution: a Gaussian dome in sigma coordinates
###############################################################
#
# A Gaussian surface on a flat bed:
#
#     b(x, y) = 0,    z_srf(x, y) = H(x, y) = G(x, y) = H₀ + A exp(-(x² + y²) / 2s²)
#
# Pagos' sigma coordinate is ζ = (z - b) / H, zero at the bed and one at the surface, so
# a point at sigma level ζ sits at physical elevation
#
#     z(x, y, ζ) = ζ G(x, y).
#
# The field being differentiated is the *physical elevation itself*, stored on the sigma
# grid as F(x, y, ζ) = ζ G(x, y). Its gradient in untransformed (x, y, z) space is known
# with no computation at all — z is its own coordinate, so ∇z = (0, 0, 1). That makes the
# reference exact: there is no discretization error hiding in the expected answer, which
# is what "compare against the result without any transform" has to mean if the
# comparison is to be worth anything.
#
# In sigma space neither term is trivial, and both carry the Gaussian:
#
#     ∂F/∂x|_ζ = ζ Gₓ        ∂F/∂ζ = G
#
# and the chain rule for ζ = z / G, i.e. ∂ζ/∂x|_z = -(ζ/G) Gₓ, gives back
#
#     ∂z/∂x|_z = ∂F/∂x|_ζ + (∂ζ/∂x|_z) ∂F/∂ζ = ζGₓ - (ζ/G) Gₓ G = 0
#     ∂z/∂z    = (1/G) ∂F/∂ζ                                     = 1
#
# The cancellation is exact discretely too, not just analytically — but only because every
# term below is evaluated at the same node. That is the reason the diagnostics are
# collocated back onto `aa` rather than left at the `acx`/z-`Vertex` nodes the operators
# naturally produce: interpolating ζ to a face and evaluating Gₓ at one leaves a residual
# of (ζ̄ - ζ_face) Gₓ, which is nonzero on a stretched sigma axis. Tests 4 and 5 would
# then need a discretization-sized tolerance instead of a machine-precision one, and would
# stop distinguishing "the transform is right" from "the transform is nearly right".

_gaussian(x, y, p)    = p.H₀ + p.A * exp(-(x^2 + y^2) / (2 * p.s^2))
_gaussian_dx(x, y, p) = -p.A * x / p.s^2 * exp(-(x^2 + y^2) / (2 * p.s^2))

@kernel inbounds = true function _fill_gaussian!(H, grid, p, O)
    I = @index(Global, NTuple)
    I = I + O
    x, y, _ = coord(grid, location(H), I...)
    H[I...] = _gaussian(x, y, p)
end

# Physical elevation on the sigma grid. Reads the depth-integrated `H` at k == 1 while the
# launch is driven by the column grid — the two-grids-per-domain pattern of Phase 2.
@kernel inbounds = true function _fill_elevation!(F, H, grid, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, k = I
    ζ = zcoord(grid, location(F), i, j, k)
    F[i, j, k] = ζ * H[i, j, 1]
end

# Both the raw sigma-space derivatives and the transformed physical ones, so the test can
# check that the cancellation in `dzdx` is a real cancellation of two large terms rather
# than 0 - 0.
@kernel inbounds = true function _sigma_gradients!(dzdx, dzdz, Fx, Fζ, F, H, grid, grid2d, O)
    I = @index(Global, NTuple)
    I = I + O
    i, j, k = I

    # Average the operators' natural outputs (x-`Vertex`, z-`Vertex`) back onto `aa`, so
    # every term of the chain rule below lives at the same point. `∂z_σ` (not Chmy's own
    # `∂z`) because this is exactly the sigma axis it is wrong on — see
    # `src/api/sigma_operators.jl`.
    dFdx = (∂x(F, grid, i, j, k) + ∂x(F, grid, i + 1, j, k)) / 2
    dFdζ = (∂z_σ(F, grid, i, j, k) + ∂z_σ(F, grid, i, j, k + 1)) / 2
    dGdx = (∂x(H, grid2d, i, j, 1) + ∂x(H, grid2d, i + 1, j, 1)) / 2

    ζ = zcoord(grid, location(F), i, j, k)
    G = H[i, j, 1]

    Fx[i, j, k] = dFdx
    Fζ[i, j, k] = dFdζ

    dzdx[i, j, k] = dFdx - (ζ / G) * dGdx * dFdζ
    dzdz[i, j, k] = dFdζ / G
end

@testset "sigma-coordinate gradients of a Gaussian dome" begin
    lx, ly, dx, dy = 1000e3, 1000e3, 25e3, 25e3
    p = (H₀ = 1000.0, A = 2000.0, s = 200e3)

    # Both a uniform sigma axis and a stretched one. The stretched case is the one that
    # matters — sigma layering exists to refine towards the bed — and it is the case where
    # a naive vertical difference silently picks up the wrong spacing.
    for transform in (LinearSigmaTransform(Float64, 12), QuadraticSigmaTransform(Float64, 12))
        layering = CorrectedVerticalLayering(Float64, transform)
        grid     = StaggeredGrid(Float64, lx, ly, dx, dy, layering)
        rt       = Runtime(grid)
        # `Runtime` carries a Launcher for the column grid only, so the depth-integrated
        # fill needs its own (pagos-roadmaps/chmy.md §3, "Two grids per domain").
        launch2d = Launcher(grid.arch, grid.grid2d)

        @testset "$(nameof(typeof(transform))), p = $(transform.exponent)" begin
            H = Field(rt.arch, grid.grid2d, Center())
            F = Field(rt.arch, rt.grid, Center())
            dzdx, dzdz, Fx, Fζ = ntuple(_ -> Field(rt.arch, rt.grid, Center()), 4)

            launch2d(rt.arch, grid.grid2d, _fill_gaussian! => (H, grid.grid2d, p))
            rt.launch(rt.arch, rt.grid, _fill_elevation! => (F, H, rt.grid))
            rt.launch(rt.arch, rt.grid,
                      _sigma_gradients! => (dzdx, dzdz, Fx, Fζ, F, H, rt.grid, grid.grid2d))

            x = collect(grid.x)
            ζ = collect(grid.z)
            G  = [_gaussian(xi, yi, p)    for xi in x, yi in collect(grid.y)]
            Gx = [_gaussian_dx(xi, yi, p) for xi in x, yi in collect(grid.y)]

            # ── The sigma-space derivatives are what was actually computed ──────────────
            #
            # ∂F/∂ζ = G exactly: F is linear in ζ, so a two-point vertical difference
            # scaled by the spacing it actually spans is exact, on the stretched axis too.
            # This is the assertion that fails if the vertical spacing is taken at the
            # wrong location — it comes out scaled by Δζ_center/Δζ_vertex.
            @test interior(Fζ) ≈ repeat(G, 1, 1, grid.nz) rtol = 1e-14

            # ∂F/∂x|_ζ = ζGₓ, second-order accurate in dx. Not the physical gradient —
            # this is the term the transform has to cancel.
            Fx_analytic = [ζ[k] * Gx[i, j] for i in axes(Gx, 1), j in axes(Gx, 2), k in eachindex(ζ)]
            @test maximum(abs, interior(Fx) - Fx_analytic) < 5e-5   # O(dx²·A/s³), dx/s = 1/8

            # ...and it is genuinely large — |ζGₓ| peaks at ≈ A/(s√e) ≈ 6.1e-3, nine
            # orders of magnitude above the residual asserted below — so the cancellation
            # that follows is a real one between two big terms, not 0 - 0.
            @test maximum(abs, interior(Fx)) > 0.9 * maximum(abs, Fx_analytic) > 1e-3

            # ── The transformed gradients match the untransformed answer ────────────────
            #
            # ∇z = (0, 0, 1) in physical coordinates, by definition of z. Machine
            # precision, not discretization precision: the discrete cancellation is exact
            # because every term is collocated at `aa`.
            @test maximum(abs, interior(dzdx)) < 1e-12
            @test interior(dzdz) ≈ ones(grid.nx, grid.ny, grid.nz) rtol = 1e-14

            # The recovered ∂z/∂z is a real division by the local thickness, not a
            # constant that happened to survive: it is 1 across a domain where G varies
            # by more than a factor of two.
            @test maximum(G) / minimum(G) > 2
        end
    end
end

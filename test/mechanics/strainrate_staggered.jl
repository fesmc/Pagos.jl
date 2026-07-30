using Pagos
using Test

include("../test_helpers/chmy.jl")

# `analytic_like` takes fun(x, y); the sigma-dependent variant is local to this file.
function analytic_like3d(f, grid, fun)
    out = similar(interior(f))
    loc = location(f)
    for k in axes(out, 3), j in axes(out, 2), i in axes(out, 1)
        x, y, ζ = coord(grid, loc, i, j, k)
        out[i, j, k] = fun(x, y, ζ)
    end
    return out
end

# The Chmy-native, C-grid strain-rate tensor and deviatoric stress. Validated against
# analytic solutions, plus an explicit cross-check against the collocated implementation
# (which is kept, untouched, alongside it) in the one regime where both discretizations are
# exact and therefore *must* agree.

@testset "strain rate and deviatoric stress (C-grid staggered)" begin
    lx, ly, dx, dy = 8.0, 8.0, 1.0, 1.0
    nz = 6
    H0 = 2.0        # thickness; the sigma scaling is ∂/∂z = (1/H) ∂/∂ζ

    function setup(T = Float64; nlayers = nz, thickness = H0)
        layering = CorrectedVerticalLayering(T, QuadraticSigmaTransform(T, nlayers))
        grid = StaggeredGrid(T, lx, ly, dx, dy, layering)
        rt   = Runtime(grid)
        mech, mat = MechanicState(grid), MaterialState(grid)
        fill_analytic!(mech.topography.thickness, rt.grid2d, (x, y) -> T(thickness))
        return grid, rt, mech, mat
    end

    momentum = BlatterPattynMomentumBalance()

    @testset "the layout closes: every gradient lands where its component lives" begin
        grid, _, mech, _ = setup()
        (; velocity, strainrate) = mech

        # ε̇_xx = ∂u/∂x, both at aa
        @test location(velocity.x_dx) === location(strainrate.xx) === Pagos.NODE_AA
        @test location(velocity.y_dy) === location(strainrate.yy) === Pagos.NODE_AA
        # ε̇_xy = (∂u/∂y + ∂v/∂x)/2 — all three at ab, no interpolation in the formula
        @test location(velocity.x_dy) === location(velocity.y_dx) ===
              location(strainrate.xy) === Pagos.NODE_AB
        # ε̇_xz = (∂u/∂z + ∂w/∂x)/2 — all three at acx/z-Vertex
        @test location(velocity.x_dz) === location(velocity.z_dx) ===
              location(strainrate.xz) === Pagos.NODE_ACX_AC
        @test location(velocity.y_dz) === location(velocity.z_dy) ===
              location(strainrate.yz) === Pagos.NODE_ACY_AC
        # ∂w/∂z closes the tensor back at aa
        @test location(velocity.z_dz) === location(strainrate.zz) === Pagos.NODE_AA

        # ...and the four node classes really are four different shapes, which is why a
        # shared flat linear index cannot work here.
        shapes = size.((strainrate.xx, strainrate.xy, strainrate.xz, strainrate.yz))
        @test length(unique(shapes)) == 4
    end

    # A linear velocity field: every derivative is exact for both the staggered two-point
    # difference and the collocated central difference, so the tensor is exact.
    @testset "exact: linear velocity field" begin
        a, b, c, d, e, f = 2.0, 3.0, -1.0, 0.5, 0.25, -0.75
        _, rt, mech, _ = setup()

        # u = a·x + b·y, v = c·x + d·y, w = (e·x + f·y)·ζ·H  ⟹  ∂w/∂z = e·x + f·y
        fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> a * x + b * y)
        fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> c * x + d * y)
        fill_analytic3d!(mech.velocity.z, rt.grid, (x, y, ζ) -> (e * x + f * y) * ζ * H0)

        raw_strainrate!(mech, momentum, rt)
        (; strainrate, velocity) = mech

        @test all(interior(velocity.x_dx) .≈ a)
        @test all(interior(velocity.x_dy) .≈ b)
        @test all(interior(velocity.y_dx) .≈ c)
        @test all(interior(velocity.y_dy) .≈ d)

        @test all(interior(strainrate.xx) .≈ a)
        @test all(interior(strainrate.yy) .≈ d)
        # BlatterPattyn is a FullColumnMomentumBalance ⟹ ε̇_zz from ∂w/∂z, not from
        # incompressibility. Here w was chosen so that ∂w/∂z = e·x + f·y.
        @test interior(strainrate.zz) ≈
              analytic_like3d(strainrate.zz, rt.grid, (x, y, ζ) -> e * x + f * y)
        @test all(interior(strainrate.xy) .≈ (b + c) / 2)

        # the symmetric duplicates are filled, not left at zero
        @test interior(strainrate.yx) == interior(strainrate.xy)
        @test interior(strainrate.zx) == interior(strainrate.xz)
        @test interior(strainrate.zy) == interior(strainrate.yz)
    end

    # The sigma scaling `∂/∂z = (1/H) ∂/∂ζ` on the *stretched* axis: this is where Chmy's
    # own `∂z` would be wrong and `∂z_σ` is required. u linear in ζ ⟹ ∂u/∂z = slope/H,
    # exactly, at every interface including the bed and surface.
    @testset "exact: vertical shear on the sigma axis (∂z_σ and the 1/H scaling)" begin
        slope = 7.0
        for thickness in (1.0, 250.0)
            _, rt, mech, _ = setup(; thickness)

            fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> slope * ζ)
            fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> -2slope * ζ)
            fill_analytic3d!(mech.velocity.z, rt.grid, (x, y, ζ) -> 0.0)

            raw_strainrate!(mech, momentum, rt)

            @test all(interior(mech.velocity.x_dz) .≈ slope / thickness)
            @test all(interior(mech.velocity.y_dz) .≈ -2slope / thickness)
            @test all(interior(mech.strainrate.xz) .≈ slope / (2 * thickness))
            @test all(interior(mech.strainrate.yz) .≈ -slope / thickness)
        end
    end

    # Chmy's raw `∂z` would give a visibly different (wrong) answer mid-column on the
    # stretched axis — pinned so this test is known to actually discriminate, rather than
    # passing because the axis happens to be uniform.
    @testset "the sigma axis really is stretched (so ∂z_σ matters)" begin
        layering = CorrectedVerticalLayering(Float64, QuadraticSigmaTransform(Float64, nz))
        spacings = diff(layering.ζ_aa)
        @test maximum(spacings) / minimum(spacings) > 2
    end

    @testset "no ice: the vertical shear is zero, not Inf" begin
        _, rt, mech, _ = setup(; thickness = 0.0)
        fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> 3ζ)
        fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> 3ζ)
        fill_analytic3d!(mech.velocity.z, rt.grid, (x, y, ζ) -> 0.0)

        raw_strainrate!(mech, momentum, rt)

        @test all(iszero, interior(mech.velocity.x_dz))
        @test !Pagos.hasnan(mech.strainrate.effective)
        @test all(isfinite, interior(mech.strainrate.effective))
    end

    # ε̇_zz from incompressibility for the balances that do not resolve w.
    @testset "ε̇_zz: incompressibility vs. resolved ∂w/∂z" begin
        a, d = 2.0, 3.0
        for mb in (SSAMomentumBalance(), DIVAMomentumBalance())
            _, rt, mech, _ = setup()
            fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> a * x)
            fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> d * y)
            fill_analytic3d!(mech.velocity.z, rt.grid, (x, y, ζ) -> 99.0 * ζ * H0)

            raw_strainrate!(mech, mb, rt)
            # ∂w/∂z would be 99 here; incompressibility must give -(a + d) instead
            @test all(interior(mech.strainrate.zz) .≈ -(a + d))
        end

        _, rt, mech, _ = setup()
        fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> a * x)
        fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> d * y)
        fill_analytic3d!(mech.velocity.z, rt.grid, (x, y, ζ) -> 99.0 * ζ * H0)
        raw_strainrate!(mech, BlatterPattynMomentumBalance(), rt)
        @test all(interior(mech.strainrate.zz) .≈ 99.0)
    end

    # The effective strain rate needs the off-diagonals interpolated back to `aa`, so it
    # cannot share the component pass. With a linear velocity field every component is
    # spatially constant, so the interpolation is exact and the invariant is analytic.
    @testset "effective strain rate: second invariant at aa" begin
        a, b, c, d = 2.0, 3.0, -1.0, 0.5
        _, rt, mech, _ = setup()

        fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> a * x + b * y)
        fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> c * x + d * y)
        fill_analytic3d!(mech.velocity.z, rt.grid, (x, y, ζ) -> 0.0)

        raw_strainrate!(mech, momentum, rt)

        exx, eyy, ezz = a, d, 0.0
        exy = (b + c) / 2
        expected = sqrt((exx^2 + eyy^2 + ezz^2) / 2 + exy^2)
        @test all(interior(mech.strainrate.effective) .≈ expected)
        @test location(mech.strainrate.effective) === Pagos.NODE_AA
    end

    @testset "deviatoric stress: τ = 2ηε̇ with a harmonically staggered viscosity" begin
        a, b, c, d, η0 = 2.0, 3.0, -1.0, 0.5, 1.0e13
        _, rt, mech, mat = setup()

        fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> a * x + b * y)
        fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> c * x + d * y)
        fill_analytic3d!(mech.velocity.z, rt.grid, (x, y, ζ) -> 0.0)
        fill_analytic3d!(mat.eta_ice, rt.grid, (x, y, ζ) -> η0)

        raw_strainrate!(mech, momentum, rt)
        deviatoric_stress!(mech, mat, rt)
        (; stress, strainrate) = mech

        # Uniform η ⟹ the harmonic mean is η itself, so every component is exactly 2ηε̇.
        @test interior(stress.xx) ≈ 2η0 .* interior(strainrate.xx)
        @test interior(stress.xy) ≈ 2η0 .* interior(strainrate.xy)
        @test interior(stress.xz) ≈ 2η0 .* interior(strainrate.xz)
        @test interior(stress.yz) ≈ 2η0 .* interior(strainrate.yz)
        # traceless identity
        @test interior(stress.zz) ≈ -(interior(stress.xx) .+ interior(stress.yy))

        expected = sqrt((interior(stress.xx)[1] ^ 2 + interior(stress.yy)[1] ^ 2 +
                         interior(stress.zz)[1] ^ 2) / 2 + interior(stress.xy)[1] ^ 2)
        @test all(interior(stress.effective) .≈ expected)
    end

    # The viscosity interpolation is harmonic, not arithmetic — a distinction that only
    # shows up with a contrast, and that matters at margins where η jumps.
    @testset "viscosity staggering is the harmonic mean, not the arithmetic one" begin
        _, rt, mech, mat = setup()
        η_lo, η_hi = 1.0e12, 1.0e14

        # a step in η across the x-faces, uniform in y and z
        fill_analytic3d!(mat.eta_ice, rt.grid,
                        (x, y, ζ) -> x < 0 ? η_lo : η_hi)
        fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> 0.0)
        fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> 3.0 * x)
        fill_analytic3d!(mech.velocity.z, rt.grid, (x, y, ζ) -> 0.0)

        raw_strainrate!(mech, momentum, rt)
        deviatoric_stress!(mech, mat, rt)

        # ε̇_xy is uniform (= 3/2), so τ_xy at ab reads off the staggered η directly.
        harmonic   = 2 / (1 / η_lo + 1 / η_hi)
        arithmetic = (η_lo + η_hi) / 2

        # The x-face that *straddles* the step: Chmy puts vertex `i` between centres `i-1`
        # and `i`, so the face to pick is the first one whose right-hand cell centre is
        # past x = 0 — not the first vertex past x = 0, whose two cells are both on the
        # high side. Asserted below rather than assumed.
        xc     = collect(xcenters(rt.grid))
        i_step = findfirst(>(0.0), xc)
        η_i    = interior(mat.eta_ice)
        @test η_i[i_step - 1, 1, 1] ≈ η_lo
        @test η_i[i_step, 1, 1] ≈ η_hi

        τ_face = interior(mech.stress.xy)[i_step, 2, 1]
        exy    = interior(mech.strainrate.xy)[i_step, 2, 1]

        @test exy ≈ 3 / 2
        @test τ_face ≈ 2 * harmonic * exy
        @test !isapprox(τ_face, 2 * arithmetic * exy; rtol = 1e-3)
        @test harmonic < arithmetic         # the whole point: no stiff-cell domination

        # ...and a face well inside the uniform region is just 2ηε̇, so the harmonic mean
        # is not quietly perturbing values away from the contrast.
        @test interior(mech.stress.xy)[i_step + 2, 2, 1] ≈ 2 * η_hi * exy
    end

    # A zero viscosity gives NaN through `hlerp`, not zero — documented on the method, and
    # pinned here so it is a known property rather than a surprise.
    @testset "a zero viscosity gives NaN through hlerp" begin
        _, rt, mech, mat = setup()      # mat.eta_ice is freshly allocated: all zeros
        fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> 3.0 * y)
        fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> 0.0)
        fill_analytic3d!(mech.velocity.z, rt.grid, (x, y, ζ) -> 0.0)

        raw_strainrate!(mech, momentum, rt)
        deviatoric_stress!(mech, mat, rt)
        @test Pagos.hasnan(mech.stress.xy)
    end

    ###########################################################
    # Old vs. new
    ###########################################################

    # The collocated implementation is kept and still works. In a *linear* velocity field
    # both discretizations are exact, so they must agree to roundoff despite evaluating at
    # different points — the strongest available cross-check between the two code paths,
    # and the only regime where "they disagree" would be a genuine bug rather than the two
    # schemes simply discretizing different things.
    @testset "old vs. new: agree exactly on a linear velocity field" begin
        a, b, c, d = 2.0, 3.0, -1.0, 0.5

        # --- new (staggered) ---
        _, rt, mech, _ = setup()
        fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> a * x + b * y)
        fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> c * x + d * y)
        fill_analytic3d!(mech.velocity.z, rt.grid, (x, y, ζ) -> 0.0)
        raw_strainrate!(mech, SSAMomentumBalance(), rt)

        # --- old (collocated, plain arrays on a RegularGrid) ---
        rgrid = RegularGrid(Float64, lx, ly, dx, dy)
        rmech = MechanicState(rgrid)
        (; nx, ny) = rgrid
        for (i, x) in enumerate(rgrid.x), (j, y) in enumerate(rgrid.y)
            rmech.velocity.x[i, j, 1] = a * x + b * y
            rmech.velocity.y[i, j, 1] = c * x + d * y
        end
        velocitygradients!(rmech.velocity.x_dx, rmech.velocity.x_dy,
                           rmech.velocity.y_dx, rmech.velocity.y_dy,
                           rmech.velocity.x, rmech.velocity.y, rgrid.dx, rgrid.dy)
        raw_strainrate!(rmech.strainrate, rmech.velocity, SSAMomentumBalance())

        # Both are exact, at different points, so compare against the shared analytic
        # values rather than element-by-element (the grids do not even have the same size —
        # that is the documented `RegularGrid`/`StaggeredGrid` off-by-one).
        old_interior(A) = A[2:(end - 1), 2:(end - 1), :]    # drop the one-sided edges
        for (comp, exact) in ((:xx, a), (:yy, d), (:xy, (b + c) / 2), (:zz, -(a + d)))
            new_vals = interior(getproperty(mech.strainrate, comp))
            old_vals = old_interior(getproperty(rmech.strainrate, comp))
            @test all(new_vals .≈ exact)
            @test all(old_vals .≈ exact)
            @test maximum(abs, new_vals .- exact) ≈ 0 atol = 1e-12
            @test maximum(abs, old_vals .- exact) ≈ 0 atol = 1e-12
        end

        # ...and the effective strain rate, which the two build differently (the staggered
        # one has to interpolate ε̇_xy back to `aa` first).
        expected = sqrt((a^2 + d^2 + (a + d)^2) / 2 + ((b + c) / 2)^2)
        @test all(interior(mech.strainrate.effective) .≈ expected)
        @test all(old_interior(rmech.strainrate.effective) .≈ expected)
    end

    @testset "old vs. new: the collocated path is untouched and still dispatches" begin
        rgrid = RegularGrid(Float64, lx, ly, dx, dy)
        rmech, rmat = MechanicState(rgrid), MaterialState(rgrid)
        fill!(rmat.eta_ice, 1.0e13)
        fill!(rmech.strainrate.xx, 1.0e-3)

        # 3-arg collocated methods still resolve, and the fused array kernel still runs
        @test raw_strainrate!(rmech.strainrate, rmech.velocity, momentum) === nothing
        @test deviatoric_stress!(rmech, rmat) === nothing
        @test velocitygradients!(rmech.velocity, rgrid.dx, rgrid.dy) === nothing
    end

    # ...while the same collocated call on a Field-based state is refused rather than
    # silently mixing locations.
    @testset "old vs. new: the collocated fused kernel is refused on Fields" begin
        _, _, mech, mat = setup()
        @test_throws ErrorException deviatoric_stress!(mech, mat)
    end

    @testset "Float32 end to end" begin
        _, rt, mech, mat = setup(Float32)
        fill_analytic3d!(mech.velocity.x, rt.grid, (x, y, ζ) -> 2.0f0 * x + 3.0f0 * y)
        fill_analytic3d!(mech.velocity.y, rt.grid, (x, y, ζ) -> -1.0f0 * x)
        fill_analytic3d!(mech.velocity.z, rt.grid, (x, y, ζ) -> 0.0f0)
        fill_analytic3d!(mat.eta_ice, rt.grid, (x, y, ζ) -> 1.0f13)

        raw_strainrate!(mech, momentum, rt)
        deviatoric_stress!(mech, mat, rt)

        @test eltype(mech.strainrate.xy) === Float32
        @test eltype(interior(mech.strainrate.effective)) === Float32
        @test eltype(interior(mech.stress.xy)) === Float32
        @test all(interior(mech.strainrate.xx) .≈ 2.0f0)
        @test all(interior(mech.strainrate.xy) .≈ (3.0f0 - 1.0f0) / 2)
        @test all(isfinite, interior(mech.stress.effective))
    end
end

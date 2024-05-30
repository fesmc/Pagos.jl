using Pagos
using Test
include("../test_helpers/structs.jl")
include("../test_helpers/utils.jl")

function slab_icesheet(slab, dx; nx = 11, ny = 3, dtau_scaling = 1.0)
    T = Float64
    lx = nx * dx
    ly = ny * dx
    dy = dx
    domain = Domain(T, lx, ly, dx, dy)
    state = State(domain)
    params = Params{T}()
    options = Options{T}(maxiter = 10_000, printout_every = 1000, dtau_scaling = dtau_scaling)
    icesheet = IceSheet(state, domain, params, options)
    
    (; state, domain, params, options) = icesheet
    state.H .= slab.H
    state.mu .= slab.mu
    state.beta .= slab.beta
    for j = 1:ny
        state.z_b[:, j] = [(10000.0 - slab.alpha * x) for x in domain.x]
    end
    return icesheet
end

######################################

@testset "slab experiments" begin
    # Case 1 #
    slab1 = TestSlab()
    icesheet1 = slab_icesheet(slab1, 5e3)
    pseudo_transient!(icesheet1)
    @test isapprox(icesheet1.state.ux,
        fill(slab1.ub, icesheet1.domain.nx, icesheet1.domain.ny), rtol = 0.001)

    # Case 2 #
    slab2 = TestSlab(H = 500.0, mu = 4e5, beta = 30.0, alpha = 1e-3)
    icesheet2 = slab_icesheet(slab2, 5e3, dtau_scaling = 10.0)
    pseudo_transient!(icesheet2)
    @test isapprox(icesheet2.state.ux,
        fill(slab2.ub, icesheet2.domain.nx, icesheet2.domain.ny), rtol = 0.001)
end
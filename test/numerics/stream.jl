using Pagos
using Test
include("../test_helpers/structs.jl")
include("../test_helpers/utils.jl")

function solve_stream(stream, dx; nx = 11, ny = 3, dtau_scaling = 1.0)

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
    state.H .= stream.H
    state.mu .= stream.mu
    state.beta .= stream.beta
    for j = 1:ny
        state.z_b[:, j] = [(10000.0 - stream.alpha * x) for x in domain.x]
    end
    pseudo_transient!(icesheet)
    return icesheet
end

######################################

@testset "stream experiments" begin
    # Case 1 #
    stream1 = TestSlab()
    sol1 = solve_stream(stream1, 5e3)
    @test isapprox(sol1.state.ux, fill(8.93, sol1.domain.nx, sol1.domain.ny), rtol = 0.01)

    # Case 2 #
    stream2 = TestSlab(H = 500.0, mu = 4e5, beta = 30.0, alpha = 1e-3)
    sol2 = solve_stream(stream2, 5e3, dtau_scaling = 10.0)
    @test isapprox(sol2.state.ux, fill(149.0, sol2.domain.nx, sol2.domain.ny), rtol = 0.01)
end
using Pagos
using BenchmarkTools
include("../test/test_helpers/structs.jl")
include("../test/test_helpers/utils.jl")
include("../test/numerics/slab.jl")

function ice_cap(domain, r, h)
    R = sqrt.(domain.X .^ 2 + domain.Y .^ 2)
    H = - 1e-3 .* R .^ 2 .+ h
    return max.(H, 0)
end

T = Float64
domain = Domain(T, 6000, 6000, 16, 16)
state = State(domain)
params = Params{T}()
options = Options{T}(maxiter = 1000)
icesheet = IceSheet(state, domain, params, options)

(; state, domain, params, options) = icesheet
state.H = ice_cap(domain, 3000, 1000)
state.f_ice .= T.(state.H .> 0)
state.mu .= 1e5
state.β .= 1e3
state.β_acx .= state.β
state.β_acy .= state.β

@btime pseudo_dotvel!($state, $domain, $params, $options)
# 4.638 ms (0 allocations: 0 bytes)

dtau = 1.0
@btime pseudo_vel!($state.ux, $state.ux_old,
    $state.dotvel_x, $dtau, $options.theta_v)
# 47.941 μs (0 allocations: 0 bytes)


slab1 = TestSlab()
sol1 = solve_slab(slab1, 500e3, nx = 350, ny = 350, dtau_scaling = 1e-3)
@btime solve_slab($slab1, 500e3, nx = 350, ny = 350, dtau_scaling = 1e-3)
# pseudo-transient: 86.182 ms (369089 allocations: 87.11 MiB)

# In comparison, the linear solver gives: 3.573 s (1066 allocations: 1010.89 MiB)
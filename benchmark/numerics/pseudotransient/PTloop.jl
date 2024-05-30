using Pagos
using BenchmarkTools
include("../../../test/test_helpers/structs.jl")
include("../../../test/test_helpers/utils.jl")
include("../../../test/numerics/slab.jl")
# include("utils.jl")

"""
    sample_resolutions(range::Tuple{Int, Int}, Lx)

Sample resolutions for a given range of resolutions (powers of 2) and domain size (in m).
"""
function sample_resolutions(range::Tuple{Int, Int}, Lx)
    @assert range[1] < range[2] "Invalid range"
    power_vec = reverse(range[1]:range[2])
    dx_vec = 2 .^ power_vec .* 1e3
    nx_vec = round.(Int, Lx ./ dx_vec)
    return power_vec, dx_vec, nx_vec
end

function btime_slab_problem(icesheet)
    (; state, domain, params, options) = icesheet
    pseudodotvel_time = @belapsed pseudo_dotvel!($state, $domain, $params, $options)
    # 4.638 ms (0 allocations: 0 bytes)

    dtau = 1.0
    pseudovel_time = @belapsed pseudo_vel!($state.ux, $state.ux_old,
        $state.dotvel_x, $dtau, $options.theta_v)
    # 47.941 μs (0 allocations: 0 bytes)

    solve_time = @belapsed pseudo_transient!($icesheet)
    # pseudo-transient: 86.182 ms (369089 allocations: 87.11 MiB)

    return pseudodotvel_time, pseudovel_time, solve_time
end

function time_slab_problem(icesheet)
    (; state, domain, params, options) = icesheet
    dtau = 1.0

    pseudodotvel_time = @elapsed pseudo_dotvel!(state, domain, params, options)
    pseudovel_time = @elapsed pseudo_vel!(state.ux, state.ux_old,
        state.dotvel_x, dtau, options.theta_v)
    solve_time = @elapsed pseudo_transient!(icesheet)
    return pseudodotvel_time, pseudovel_time, solve_time
end
include(joinpath( @__DIR__, "../src/Pagos.jl"))
using .Pagos

T = Float64
domain = Domain(T, 6000, 6000, 16, 16)
sinx = lazy_map(sin, domain.X)
cossinx = lazy_map(cos, sinx)
ph = copy(domain.null)
@btime ph .= cossinx .+ 1.0
# 1.814 ms (1 allocation: 64 bytes)

@btime ph .= cos.(sin.(domain.X)) .+ 1.0
# 1.763 ms (4 allocations: 160 bytes)
include(joinpath( @__DIR__, "../src/Pagos.jl"))
using .Pagos
# using CairoMakie

function ice_cap(domain, r, h)
    R = sqrt.(domain.X .^ 2 + domain.Y .^ 2)
    H = - 1e-3 .* R .^ 2 .+ h
    return max.(H, 0)
end

T = Float64
domain = Domain(T, 6000, 6000, 16, 16)

state = State(domain)
params = Params{T}()
options = Options{T}()
icesheet = IceSheet(state, domain, params, options)

(; state, domain, params, options) = icesheet
state.H = ice_cap(domain, 3000, 1000)
# fig, ax, srf = surface(domain.X, domain.Y, state.H)

state.f_ice .= T.(state.H .> 0)
state.mu .= 1e5
state.β .= 1e3
state.β_acx .= state.β
state.β_acy .= state.β

pseudo_dotvel!(state, domain, params)
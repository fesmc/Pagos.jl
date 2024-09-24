"""
    advect!(icesheet::IceSheet)

Update the ice thickness field using the advection equation.
"""
function advect!(icesheet::IceSheet)
    return advect!(icesheet.state, icesheet.domain, icesheet.options)
end

# TODO: handle the new ice-covered cells (otherwise, PT crashes)
function advect!(state::State{T}, domain::Domain{T}, options::Options{T}) where {T<:AbstractFloat}
    (; H, ux, uy, ux_old, uy_old, prealloc) = state
    (; dx, dy, nx, ny) = domain
    dt = 1.0
    delx!(prealloc, H .* ux, dx, nx)
    @. H += prealloc * dt
    dely!(prealloc, H .* uy, dy, ny)
    @. H += prealloc * dt

    ux[H .< 0] .= 0.0
    uy[H .< 0] .= 0.0
    ux_old .= ux
    uy_old .= uy
    H[H .< 0] .= 0.0
    return nothing
end
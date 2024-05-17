mutable struct Tools{T<:AbstractFloat}
    n_terms::Int
    n_u::Int
    n_sprs::Int
    u::Vector{T} = fill(0.0, n_u)
    b::Vector{T} = fill(0.0, n_u)
    Ai::Vector{Int} = fill(0, n_sprs)
    Aj::Vector{Int} = fill(0, n_sprs)
    Av::Vector{T} = fill(0.0, n_sprs)
end

function Tools(nx::Int, ny::Int; n_terms::Int = 9)
    n_u = 2 * nx * ny
    n_sprs = n_u * n_terms

    # Populate dense array for now
    u = fill(0.0, n_u)
    b = fill(0.0, n_u)

    # Define array A component vectors (I, J, V)
    Ai = fill(0, n_sprs)
    Aj = fill(0, n_sprs)
    Av = fill(0.0, n_sprs)
    return Tools(n_terms, n_u, n_sprs, u, b, Ai, Aj, Av)
end

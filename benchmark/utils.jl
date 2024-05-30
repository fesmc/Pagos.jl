function ice_cap(domain, r, h)
    R = sqrt.(domain.X .^ 2 + domain.Y .^ 2)
    H = - 1e-3 .* R .^ 2 .+ h
    return max.(H, 0)
end
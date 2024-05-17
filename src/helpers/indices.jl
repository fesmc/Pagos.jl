function periodic_bc_minusindex(i::Int, n::Int)
    im1 = i - 1
    if im1 == 0
        im1 = n
    end
    return im1
end

function periodic_bc_plusindex(i::Int, n::Int)
    ip1 = i + 1
    if ip1 == n + 1
        ip1 = 1
    end
    return ip1
end

periodic_bc_indices(i, n) = periodic_bc_minusindex(i, n)..., periodic_bc_plusindex(i, n)...
periodic_bc_indices(ix, iy, nx, ny) =
    periodic_bc_indices(ix, nx)..., periodic_bc_indices(iy, ny)...

function ij2n_ux(i, j, nx, ny)
    n = (i - 1) * ny + j
    return n
end

function ij2n_uy(i, j, nx, ny)
    n = (i - 1) * ny + j + nx * ny
    return n
end

function halfbound(x, i::Int, im1::Int, ip1::Int)
    if x > 0
        i0 = i
        i1 = ip1
        wt = x / 2.0
    else
        i0 = im1
        i1 = i
        wt = 1.0 - abs(x) / 2.0
    end
    return i0, i1, wt
end
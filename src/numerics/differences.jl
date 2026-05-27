"""
    delx1(u, dx)
    delx1(u, dx, nx)

Compute the derivative of `u` with respect to `x` using a central difference scheme.
"""
function delx1(u, dx)
    nx = size(u, 1)
    return delx1(u, dx, nx)
end

function delx1(u, dx, mask)
    du = similar(u)
    delx1!(du, u, dx, mask)
    return du
end

function delx1!(du, u, dx, mask)
    @inbounds for I in CartesianIndices(mask)[mask]
        i, j = Tuple(I)
        if mask[i+1, j] && mask[i-1, j]
            du[i, j] = (u[i+1, j] - u[i-1, j]) / (dx * 2)
        elseif mask[i+1, j]
            du[i, j] = (u[i+1, j] - u[i, j]) / dx
        elseif mask[i-1, j]
            du[i, j] = (u[i, j] - u[i-1, j]) / dx
        else
            du[i, j] = 0.0
        end
    end
end

"""
    delx1!(du, u, dx, nx)

Update `du` (in place), the derivative of `u` with respect to `x` using a central
difference scheme.
"""
# function delx1!(du, u, dx, nx::Int)
#     @inbounds for j in axes(du, 2)
#         for i in axes(du, 1)[2:nx-1]
#             du[i, j] = (u[i+1, j] - u[i-1, j]) / (dx * 2)
#         end
#         du[1, j] = (u[2, j] - u[1, j]) / dx
#         du[nx, j] = (u[nx, j] - u[nx-1, j]) / dx
#     end
# end

"""
    delx2(u, dy)
    delx2(u, dy, ny)

Compute the derivative of `u` with respect to `y` using a central difference scheme.
"""
function delx2(u, dy)
    ny = size(u, 2)
    return delx2(u, dy, ny)
end

function delx2(u, dy, mask)
    du = similar(u)
    delx2!(du, u, dy, mask)
    return du
end

function delx2!(du, u, dy, mask)
    @inbounds for I in CartesianIndices(mask)[mask]
        i, j = Tuple(I)
        if mask[i, j+1] && mask[i, j-1]
            du[i, j] = (u[i, j+1] - u[i, j-1]) / (dy * 2)
        elseif mask[i, j+1]
            du[i, j] = (u[i, j+1] - u[i, j]) / dy
        elseif mask[i, j-1]
            du[i, j] = (u[i, j] - u[i, j-1]) / dy
        else
            du[i, j] = 0.0
        end
    end
end

"""
    delx2!(du, u, dy, ny)

Update `du` (in place), the derivative of `u` with respect to `y` using a central
difference scheme.
"""
# function delx2!(du, u, dy, ny::Int)
#     @inbounds for i in axes(du, 1)
#         for j in axes(du, 2)[2:ny-1]
#             du[i, j] = (u[i, j+1] - u[i, j-1]) / (dy * 2)
#         end
#         du[i, 1] = (u[i, 2] - u[i, 1]) / dy
#         du[i, ny] = (u[i, ny] - u[i, ny-1]) / dy
#     end
# end

"""
Derivatives accounting for the vertical coordinate transformation.
"""
# TODO: the four code lines below are only improved pseudo-code. Needs
# to be well implemented!
# Commented out to avoid CI test failure (function signatures below conflict with those above)
#delx1(u, delzt_delx1) = delx1t(u) + delzt_delx1 * delzt(u)
#delx2(u, delzt_delx2) = delx2t(u) + delzt_delx2 * delzt(u)
#delz(u, delzt_delz) = delzt_delz * delzt(u)
#delt(u, delzt_delt) = deltt(u) + delzt_delt * delzt(u)
"""
    delx(u, dx)
    delx(u, dx, nx)

Compute the derivative of `u` with respect to `x` using a central difference scheme.
"""
function delx(u, dx)
    nx = size(u, 1)
    return delx(u, dx, nx)
end

function delx(u, dx, nx)
    du = similar(u)
    delx!(du, u, dx, nx)
    return du
end

"""
    delx!(du, u, dx, nx)

Update `du` (in place), the derivative of `u` with respect to `x` using a central
difference scheme.
"""
function delx!(du::Matrix{T}, u::Matrix{T}, dx::T, nx::Int) where {T<:AbstractFloat}
    @inbounds for j in axes(du, 2)
        for i in axes(du, 1)[2:nx-1]
            du[i, j] = (u[i+1, j] - u[i-1, j]) / (dx * 2)
        end
        du[1, j] = (u[2, j] - u[1, j]) / dx
        du[nx, j] = (u[nx, j] - u[nx-1, j]) / dx
    end
end

"""
    dely(u, dy)
    dely(u, dy, ny)

Compute the derivative of `u` with respect to `y` using a central difference scheme.
"""
function dely(u, dy)
    ny = size(u, 2)
    return dely(u, dy, ny)
end

function dely(u, dy, ny)
    du = similar(u)
    dely!(du, u, dy, ny)
    return du
end

"""
    dely!(du, u, dy, ny)

Update `du` (in place), the derivative of `u` with respect to `y` using a central
difference scheme.
"""
function dely!(du::Matrix{T}, u::Matrix{T}, dy::T, ny::Int) where {T<:AbstractFloat}
    @inbounds for i in axes(du, 1)
        for j in axes(du, 2)[2:ny-1]
            du[i, j] = (u[i, j+1] - u[i, j-1]) / (dy * 2)
        end
        du[i, 1] = (u[i, 2] - u[i, 1]) / dy
        du[i, ny] = (u[i, ny] - u[i, ny-1]) / dy
    end
end
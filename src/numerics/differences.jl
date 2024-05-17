function delx(u, dx)
    du = similar(u)
    nx = size(u, 1)
    delx!(du, u, dx, nx)
    return du
end

function delx!(du::Matrix{T}, u::Matrix{T}, dx::T, nx::Int) where {T<:AbstractFloat}
    @inbounds for j in axes(du, 2)
        for i in axes(du, 1)[2:nx-1]
            du[i, j] = (u[i+1, j] - u[i-1, j]) / (dx * 2)
        end
        du[1, j] = (u[2, j] - u[1, j]) / dx
        du[nx, j] = (u[nx, j] - u[nx-1, j]) / dx
    end
end

function dely(u, dy)
    du = similar(u)
    ny = size(u, 2)
    dely!(du, u, dy, ny)
    return du
end

function dely!(du::Matrix{T}, u::Matrix{T}, dy::T, ny::Int) where {T<:AbstractFloat}
    @inbounds for i in axes(du, 1)
        for j in axes(du, 2)[2:ny-1]
            du[i, j] = (u[i, j+1] - u[i, j-1]) / (dy * 2)
        end
        du[i, 1] = (u[i, 2] - u[i, 1]) / dy
        du[i, ny] = (u[i, ny] - u[i, ny-1]) / dy
    end
end
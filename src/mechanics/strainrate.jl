"""
    scaledstrainrate!(strainrate_xx, strainrate_xy, strainrate_yy, ux_x, ux_y, uy_x, uy_y, N_ab)

Compute the strain rate tensor components scaled by `` 2 \\, \\mu \\, H ``.
"""
function scaledstrainrate!(strainrate_xx::M, strainrate_xy::M, strainrate_yy::M,
    ux_x::M, ux_y::M, uy_x::M, uy_y::M, N_ab::M) where {T<:AbstractFloat, M<:Matrix{T}}
    @. strainrate_xx = T(2.0) * N_ab * (T(2.0) * ux_x + uy_y)
    @. strainrate_xy = N_ab * (ux_y + uy_x)
    @. strainrate_yy = T(2.0) * N_ab * (ux_x + T(2.0) * uy_y)
    return nothing
end

# function strainrate!(strainrate_xx, strainrate_xy, strainrate_yy, ux_x, ux_y, uy_x, uy_y)
#     @. strainrate_xx = 2.0 * ux_x + uy_y
#     @. strainrate_xy = 0.5 .* (ux_y + uy_x)
#     @. strainrate_yy = ux_x + 2.0 * uy_y
#     return nothing
# end

"""
    velocitygradients!(ux_x, ux_y, uy_x, uy_y, ux, uy, dx, dy, nx, ny)

Compute the velocity gradients in x (`ux_x, uy_x`) and y-direction (`ux_y, uy_y`).
The gradients are computed using the central difference scheme. The input velocities
`ux` and `uy` are defined on a staggered grid with dimensions `nx` and `ny`.
The grid spacing in x and y-direction is given by `dx` and `dy`.
"""
function velocitygradients!(ux_x, ux_y, uy_x, uy_y, ux, uy, dx, dy, nx, ny)
    delx!(ux_x, ux, dx, nx)
    dely!(ux_y, ux, dy, ny)
    delx!(uy_x, uy, dx, nx)
    dely!(uy_y, uy, dy, ny)
    return nothing
end
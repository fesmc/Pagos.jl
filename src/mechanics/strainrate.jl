function strainrate!(strainrate_xx, strainrate_xy, strainrate_yy, ux_x, ux_y, uy_x, uy_y, mu, H)
    @. strainrate_xx = 2.0 * mu * H * (2.0 * ux_x + uy_y)
    @. strainrate_xy = mu * H * (ux_y + uy_x)
    @. strainrate_yy = 2.0 * mu * H * (ux_x + 2.0 * uy_y)
    return nothing
end

function velocitygradients!(ux_x, ux_y, uy_x, uy_y, ux, uy, dx, dy, nx, ny)
    delx!(ux_x, ux, dx, nx)
    dely!(ux_y, ux, dy, ny)
    delx!(uy_x, uy, dx, nx)
    dely!(uy_y, uy, dy, ny)
    return nothing
end
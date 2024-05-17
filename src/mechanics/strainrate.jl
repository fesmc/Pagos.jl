function strainrate!(strainrate_xx, strainrate_xy, strainrate_yy, ux_x, ux_y, uy_x, uy_y, mu, H)
    muH = 2 .* mu .* H
    strainrate_xx .= muH .* (2 .* ux_x + uy_y)
    strainrate_xy .= muH .* 0.5 .* (ux_y + uy_x)
    strainrate_yy .= muH .* (ux_x + 2 .* uy_y)
    return nothing
end

function velocitygradients!(ux_x, ux_y, uy_x, uy_y, ux, uy, dx, dy, nx, ny)
    delx!(ux_x, ux, dx, nx)
    dely!(ux_y, ux, dy, ny)
    delx!(uy_x, uy, dx, nx)
    dely!(uy_y, uy, dy, ny)
    return nothing
end
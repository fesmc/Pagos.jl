
"""
    acx_to_nodes(varx, i, j, xn, yn)

Interpolate a variable defined on acx-nodes to the node locations of interest.
Nodes are defined within the aa-cell where the center of the cell is (0,0), 
left border is (-1,0), right border is (1,0), etc. i, j indices refer to the current
cell, and the acx-node has (i,j) defined as the right-border node. 
"""
function acx_to_nodes(varx, i, j, xn, yn)

    n = length(xn)
    nx, ny = size(varx)
    im1, ip1, jm1, jp1 = periodic_bc_indices(i, j, nx, ny)

    # Initialize interpolated variable at nodes of interest
    varn = fill(0.0, n)

    for k = 1:n

        j0, j1, wt = halfbound(yn[k], j, jm1, jp1)

        # Get left and right-side 
        v0 = (1.0 - wt) * varx[im1, j0] + wt * varx[im1, j1]
        v1 = (1.0 - wt) * varx[i, j0] + wt * varx[i, j1]

        # Interpolate horizontally to the node location 
        wt = (1.0 + xn[k]) / 2.0
        varn[k] = (1.0 - wt) * v0 + wt * v1

    end

    return varn
end

"""
    acy_to_nodes(vary, i, j, xn, yn)

Interpolate a variable defined on acy-nodes to the node locations of interest.
Nodes are defined within the aa-cell where the center of the cell is (0,0),
left border is (0,-1), right border is (0,1), etc. i, j indices refer to the current
cell, and the acy-node has (i,j) defined as the top-border node.
"""
function acy_to_nodes(vary, i, j, xn, yn)

    n = length(xn)
    nx, ny = size(vary)
    im1, ip1, jm1, jp1 = periodic_bc_indices(i, j, nx, ny)

    # Initialize interpolated variable at nodes of interest
    varn = fill(0.0, n)

    for k = 1:n
        i0, i1, wt = halfbound(xn[k], i, im1, ip1)

        # Get left and right-side 
        v0 = (1.0 - wt) * vary[i0, jm1] + wt * vary[i1, jm1]
        v1 = (1.0 - wt) * vary[i0, j] + wt * vary[i1, j]

        # Interpolate vertically to the node location 
        wt = (1.0 + yn[k]) / 2.0
        varn[k] = (1.0 - wt) * v0 + wt * v1
    end

    return varn
end

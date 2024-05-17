function calc_visc_eff_2D_aa(
    ux,
    uy,
    ATT,
    f_ice,
    dx,
    dy;
    n_glen = 3,
    eps_0 = 1e-6,
)
    # Calculate 3D effective viscosity following L19, Eq. 2
    # Use of eps_0 ensures non-zero positive viscosity value everywhere 
    # Note: viscosity is first calculated on ab-nodes, then 
    # unstaggered back to aa-nodes. This ensures more stability for 
    # visc_eff (less likely to blow up for low strain rates). 

    visc_min = 1e5

    nx, ny = size(ux)

    # Calculate exponents 
    p1 = (1.0 - n_glen) / (2.0 * n_glen)
    p2 = -1.0 / n_glen

    # Calculate squared minimum strain rate 
    eps_0_sq = eps_0 * eps_0

    # Calculate visc_eff on aa-nodes
    visc = fill(visc_min, nx, ny)
    eps_aa = fill(eps_0_sq, nx, ny)

    for i = 1:nx
        for j = 1:ny

            if f_ice[i, j] == 1.0

                im1, ip1, jm1, jp1 = periodic_bc_indices(i, j, nx, ny)

                # Get strain rate terms
                dudx_aa = (ux[i, j] - ux[im1, j]) / dx
                dvdy_aa = (uy[i, j] - uy[i, jm1]) / dy

                dudy_aa_1 = (ux[i, jp1] - ux[i, jm1]) / (2.0 * dy)
                dudy_aa_2 = (ux[im1, jp1] - ux[im1, jm1]) / (2.0 * dy)
                dudy_aa = 0.5 * (dudy_aa_1 + dudy_aa_2)

                dvdx_aa_1 = (uy[ip1, j] - uy[im1, j]) / (2.0 * dx)
                dvdx_aa_2 = (uy[ip1, jm1] - uy[im1, jm1]) / (2.0 * dx)
                dvdx_aa = 0.5 * (dvdx_aa_1 + dvdx_aa_2)

                # Calculate the total effective strain rate from L19, Eq. 21 
                eps_sq_aa =
                    dudx_aa^2 +
                    dvdy_aa^2 +
                    dudx_aa * dvdy_aa +
                    0.25 * (dudy_aa + dvdx_aa)^2 +
                    eps_0_sq
                eps_aa[i, j] = sqrt(eps_sq_aa)

                # Get rate factor on central node
                ATT_aa = ATT[i, j]

                # Calculate effective viscosity on ab-nodes
                visc[i, j] = 0.5 * (eps_sq_aa)^(p1) * ATT_aa^(p2)

            end
        end
    end

    #println("eps: ", extrema(eps_aa))

    return visc

end

"""
    calc_visc_eff_2D_nodes(ux,uy,ATT,H_ice,f_ice,dx,dy,xn,yn;n_glen=3,eps_0=1e-6,wtn=fill(1.0,length(xn)))

Calculate 3D effective viscosity following L19, Eq. 2
Use of eps_0 ensures non-zero positive viscosity value everywhere 
Note: viscosity is first calculated on ab-nodes, then 
unstaggered back to aa-nodes. This ensures more stability for 
visc_eff (less likely to blow up for low strain rates). 

Given ux on acx-nodes and uy on acy-nodes, get both quantities 
on node locations of choice [xn;yn]. Viscosity will be calculated
at those locations and the desired weighting wtn will be applied to each node.
"""
function calc_visc_eff_2D_nodes(
    ux,
    uy,
    ATT,
    f_ice,
    dx,
    dy;
    n_glen = 3,
    eps_0 = 1e-6,
)

    visc_min = 1e5
    nx, ny = size(ux)

    # Calculate exponents 
    p1 = (1.0 - n_glen) / (2.0 * n_glen)
    p2 = -1.0 / n_glen

    # Calculate squared minimum strain rate 
    eps_0_sq = eps_0 * eps_0

    # Populate strain rates over the whole domain on acx- and acy-nodes
    dudx = fill(0.0, nx, ny)
    dvdy = fill(0.0, nx, ny)
    dudy = fill(0.0, nx, ny)
    dvdx = fill(0.0, nx, ny)

    for i = 1:nx
        for j = 1:ny
            im1, ip1, jm1, jp1 = periodic_bc_indices(i, j, nx, ny)
            dudx[i, j] = (ux[ip1, j] - ux[im1, j]) / (2.0 * dx)
            dudy[i, j] = (ux[i, jp1] - ux[i, jm1]) / (2.0 * dy)
            dvdx[i, j] = (uy[ip1, j] - uy[im1, j]) / (2.0 * dx)
            dvdy[i, j] = (uy[i, jp1] - uy[i, jm1]) / (2.0 * dy)
        end
    end

    # Calculate visc_eff on aa-nodes

    visc = fill(visc_min, nx, ny)
    eps_aa = fill(eps_0_sq, nx, ny)

    wt0 = 1.0 / sqrt(3)
    xn = [wt0, -wt0, -wt0, wt0]
    yn = [wt0, wt0, -wt0, -wt0]
    wtn = [1.0, 1.0, 1.0, 1.0]

    for i = 1:nx
        for j = 1:ny

            if f_ice[i, j] == 1.0

                im1, ip1, jm1, jp1 = periodic_bc_indices(i, j, nx, ny)

                # Get strain rate terms on node locations
                dudxn = acx_to_nodes(dudx, i, j, xn, yn)
                dudyn = acx_to_nodes(dudy, i, j, xn, yn)

                dvdxn = acy_to_nodes(dvdx, i, j, xn, yn)
                dvdyn = acy_to_nodes(dvdy, i, j, xn, yn)

                # Calculate the total effective strain rate from L19, Eq. 21 
                eps_sq_n =
                    dudxn .^ 2 + dvdyn .^ 2 .+ dudxn .* dvdyn .+
                    0.25 .* (dudyn .+ dvdxn) .^ 2 .+ eps_0_sq
                eps_aa = sum(sqrt.(eps_sq_n)) / length(eps_sq_n)

                # Get rate factor on central node
                ATT_aa = ATT[i, j]

                # Calculate effective viscosity on ab-nodes
                viscn = 0.5 .* (eps_sq_n) .^ (p1) .* ATT_aa^(p2)
                visc[i, j] = sum(viscn .* wtn) / sum(wtn)

            end
        end
    end
    #println("eps: ", extrema(eps_aa))
    return visc
end
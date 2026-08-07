"""
    acx_to_nodes(varx,i,j,xn,yn)

    Interpolate a variable defined on acx-nodes
    to the node locations of interest. Nodes are 
    defined within the aa-cell where the center of the 
    cell is (0,0), left border is (-1,0), right border is (1,0), etc. 
    i, j indices refer to the current cell, and the acx-node has 
    (i,j) defined as the right-border node. 

"""
function acx_to_nodes(varx,i,j,xn,yn)

    n = length(xn);
    
    nx, ny = size(varx);

    # BC: Periodic boundary conditions
    im1 = i-1
    if im1 == 0
        im1 = nx
    end
    ip1 = i+1
    if ip1 == nx+1
        ip1 = 1
    end

    jm1 = j-1
    if jm1 == 0
        jm1 = ny
    end
    jp1 = j+1
    if jp1 == ny+1
        jp1 = 1
    end

    # Initialize interpolated variable at nodes of interest
    varn = fill(0.0,n);

    for k = 1:n

        if yn[k] > 0
            j0 = j
            j1 = jp1
            wt = yn[k] / 2.0 
        else
            j0 = jm1
            j1 = j
            wt = 1.0 - abs(yn[k]) / 2.0
        end

        # Get left and right-side 
        v0 = (1.0-wt)*varx[im1,j0] + wt*varx[im1,j1];
        v1 = (1.0-wt)*varx[i,j0]   + wt*varx[i,j1];
        
        # Interpolate horizontally to the node location 
        wt = (1.0 + xn[k]) / 2.0;
        varn[k] = (1.0-wt)*v0 + wt*v1;

    end
    
    return varn 
end

function acy_to_nodes(vary,i,j,xn,yn)

    n = length(xn);
    
    nx, ny = size(vary);
    
    # BC: Periodic boundary conditions
    im1 = i-1
    if im1 == 0
        im1 = nx
    end
    ip1 = i+1
    if ip1 == nx+1
        ip1 = 1
    end

    jm1 = j-1
    if jm1 == 0
        jm1 = ny
    end
    jp1 = j+1
    if jp1 == ny+1
        jp1 = 1
    end

    # Initialize interpolated variable at nodes of interest
    varn = fill(0.0,n);

    for k = 1:n

        if xn[k] > 0
            i0 = i
            i1 = ip1
            wt = xn[k] / 2.0 
        else
            i0 = im1
            i1 = i
            wt = 1.0 - abs(xn[k]) / 2.0
        end

        # Get left and right-side 
        v0 = (1.0-wt)*vary[i0,jm1] + wt*vary[i1,jm1];
        v1 = (1.0-wt)*vary[i0,j]   + wt*vary[i1,j];
        
        # Interpolate vertically to the node location 
        wt = (1.0 + yn[k]) / 2.0;
        varn[k] = (1.0-wt)*v0 + wt*v1;

    end
    
    return varn 
end

# Driving stress
function calc_driving_stress(H,z_srf,dx,dy,ρ,g)
    
    nx, ny = size(H);

    taud_acx = fill(0.0,nx,ny);
    taud_acy = fill(0.0,nx,ny);

    for i in 1:nx
        for j in 1:ny

            # BC: Periodic boundary conditions in x and y 
            ip1 = i+1
            if ip1 == nx+1
                ip1 = 1 
            end
            jp1 = j+1
            if jp1 == ny+1
                jp1 = 1 
            end
            
            H_mid = 0.5*(H[i,j]+H[ip1,j]);
            taud_acx[i,j] = ρ * g * H_mid * (z_srf[ip1,j]-z_srf[i,j])/dx;

            H_mid = 0.5*(H[i,j]+H[i,jp1]);
            taud_acy[i,j] = ρ * g * H_mid * (z_srf[i,jp1]-z_srf[i,j])/dy;

        end
    end

    # BC: periodic...
    # Hack to ensure driving stress is ok in last grid point,
    # since z_srf might not be periodic...
    taud_acx[end,:] = taud_acx[end-1,:];

    return taud_acx, taud_acy
end

function stagger_beta(β)

    nx, ny = size(β);

    # Stagger β to acx and acy nodes 
    β_acx = copy(β);
    β_acy = copy(β);

    for i in 1:nx
        for j in 1:ny

            # BC: Periodic boundary conditions
            ip1 = i+1
            if ip1 == nx+1
                ip1 = 1
            end

            jp1 = j+1
            if jp1 == ny+1
                jp1 = 1
            end

            β_acx[i,j] = 0.5*(β[i,j]+β[ip1,j]);
            β_acy[i,j] = 0.5*(β[i,j]+β[i,jp1]);
        end
    end

    return β_acx, β_acy
end

"""
    calc_beta_aa_power_plastic(ux_b,uy_b,c_bed,f_ice,q,u_0)

    Calculate basal friction coefficient (beta) that
    enters the SSA solver as a function of basal velocity
    using a power-law form following Bueler and van Pelt (2015)
    
"""
function calc_beta_aa_power_plastic(ux_b::Array{T,2},uy_b::Array{T,2},
                                    c_bed::Array{T,2},f_ice::Array{T,2},q,u_0) where T
    

    # Local variables
    ub_min    = 1e-3;               # [m/yr] Minimum velocity is positive small value to avoid divide by zero
    ub_sq_min = ub_min^2;

    nx, ny = size(ux_b);

    # Initially set friction to zero everywhere
    beta = fill(0.0,nx,ny);
    
    for i = 1:nx
        for j = 1:ny

            # BC: Periodic boundary conditions
            im1 = i-1
            if im1 == 0
                im1 = nx
            end
            ip1 = i+1
            if ip1 == nx+1
                ip1 = 1
            end

            jm1 = j-1
            if jm1 == 0
                jm1 = ny
            end
            jp1 = j+1
            if jp1 == ny+1
                jp1 = 1
            end

            if f_ice[i,j] == 1.0 
                # Fully ice-covered point with some fully ice-covered neighbors 

                cb_aa = c_bed[i,j];

                if q == 1.0 
                    # Linear law, no f(ub) term

                    beta[i,j] = cb_aa / u_0;

                else
                    # Non-linear law with f(ub) term 

                    # Unstagger velocity components to aa-nodes 
                    ux_aa = 0.5*(ux_b[i,j]+ux_b[im1,j]);
                    uy_aa = 0.5*(uy_b[i,j]+uy_b[i,jm1]);
                    
                    uxy_aa = sqrt(ux_aa^2 + uy_aa^2 + ub_sq_min);
                    
                    if q == 0
                        # Plastic law

                        beta[i,j] = cb_aa * (1.0 / uxy_aa);
                        
                    else

                        beta[i,j] = cb_aa * (uxy_aa / u_0)^q * (1.0 / uxy_aa);
                    
                    end
                    
                end

            else
                # Assign minimum velocity value, no staggering for simplicity

                if q == 1.0 
                    # Linear law, no f(ub) term

                    beta[i,j] = c_bed[i,j] / u_0;

                else
                    
                    uxy_b  = ub_min;

                    if q == 0.0 
                        # Plastic law

                        beta[i,j] = c_bed[i,j] * (1.0 / uxy_b);

                    else

                        beta[i,j] = c_bed[i,j] * (uxy_b / u_0)^q * (1.0 / uxy_b);

                    end 

                end

            end

        end
    end

    return beta
    
end

function calc_beta_aa_power_plastic_nodes(ux_b::Array{T,2},uy_b::Array{T,2},
                                    c_bed::Array{T,2},f_ice::Array{T,2},q,u_0) where T
    

    # Local variables
    ub_min    = 1e-3;               # [m/yr] Minimum velocity is positive small value to avoid divide by zero
    ub_sq_min = ub_min^2;

    nx, ny = size(ux_b);

    # Initially set friction to zero everywhere
    beta = fill(0.0,nx,ny);
    
    wt0 = 1.0/sqrt(3);
    xn  = [wt0,-wt0,-wt0, wt0];
    yn  = [wt0, wt0,-wt0,-wt0];
    wtn = [1.0,1.0,1.0,1.0];

    for i = 1:nx
        for j = 1:ny

            # BC: Periodic boundary conditions
            im1 = i-1
            if im1 == 0
                im1 = nx
            end
            ip1 = i+1
            if ip1 == nx+1
                ip1 = 1
            end

            jm1 = j-1
            if jm1 == 0
                jm1 = ny
            end
            jp1 = j+1
            if jp1 == ny+1
                jp1 = 1
            end

            if f_ice[i,j] == 1.0 
                # Fully ice-covered point with some fully ice-covered neighbors 

                cb_aa = c_bed[i,j];

                if q == 1.0 
                    # Linear law, no f(ub) term

                    beta[i,j] = cb_aa / u_0;

                else
                    # Non-linear law with f(ub) term 

                    uxn  = acx_to_nodes(ux_b,i,j,xn,yn);
                    uyn  = acy_to_nodes(uy_b,i,j,xn,yn);
                    uxyn = sqrt.(uxn.^2 .+ uyn.^2 .+ ub_sq_min);

                    if q == 0
                        # Plastic law

                        betan = cb_aa .* (1.0 ./ uxyn);
                        
                    else

                        betan = cb_aa .* (uxyn ./ u_0).^q .* (1.0 ./ uxyn);
                    
                    end
                    
                    beta[i,j] = sum(betan.*wtn)/sum(wtn);

                end

            else
                # Assign minimum velocity value, no staggering for simplicity

                if q == 1.0 
                    # Linear law, no f(ub) term

                    beta[i,j] = c_bed[i,j] / u_0;

                else
                    
                    uxy_b  = ub_min;

                    if q == 0.0 
                        # Plastic law

                        beta[i,j] = c_bed[i,j] * (1.0 / uxy_b);

                    else

                        beta[i,j] = c_bed[i,j] * (uxy_b / u_0)^q * (1.0 / uxy_b);

                    end 

                end

            end

        end
    end

    return beta
    
end

function calc_visc_eff_2D_aa(ux,uy,ATT,H_ice,f_ice,dx,dy;n_glen=3,eps_0=1e-6)
    # Calculate 3D effective viscosity following L19, Eq. 2
    # Use of eps_0 ensures non-zero positive viscosity value everywhere 
    # Note: viscosity is first calculated on ab-nodes, then 
    # unstaggered back to aa-nodes. This ensures more stability for 
    # visc_eff (less likely to blow up for low strain rates). 

    visc_min = 1e5;

    nx, ny = size(ux);

    # Calculate exponents 
    p1 = (1.0 - n_glen)/(2.0*n_glen);
    p2 = -1.0/n_glen;

    # Calculate squared minimum strain rate 
    eps_0_sq = eps_0*eps_0;

    # Calculate visc_eff on aa-nodes

    visc   = fill(visc_min,nx,ny);
    eps_aa = fill(eps_0_sq,nx,ny);

    for i = 1:nx
        for j = 1:ny  

            if f_ice[i,j] == 1.0

                # BC: Periodic boundary conditions
                im1 = i-1
                if im1 == 0
                    im1 = nx
                end
                ip1 = i+1
                if ip1 == nx+1
                    ip1 = 1
                end

                jm1 = j-1
                if jm1 == 0
                    jm1 = ny
                end
                jp1 = j+1
                if jp1 == ny+1
                    jp1 = 1
                end

                # Get strain rate terms
                dudx_aa = (ux[i,j]-ux[im1,j])/dx 
                dvdy_aa = (uy[i,j]-uy[i,jm1])/dy 
                
                dudy_aa_1 = (ux[i,jp1]-ux[i,jm1])/(2.0*dy)
                dudy_aa_2 = (ux[im1,jp1]-ux[im1,jm1])/(2.0*dy)
                dudy_aa   = 0.5*(dudy_aa_1+dudy_aa_2)

                dvdx_aa_1 = (uy[ip1,j]-uy[im1,j])/(2.0*dx)
                dvdx_aa_2 = (uy[ip1,jm1]-uy[im1,jm1])/(2.0*dx)
                dvdx_aa   = 0.5*(dvdx_aa_1+dvdx_aa_2)

                # Calculate the total effective strain rate from L19, Eq. 21 
                eps_sq_aa = dudx_aa^2 + dvdy_aa^2 + dudx_aa*dvdy_aa + 0.25*(dudy_aa+dvdx_aa)^2 + eps_0_sq
                eps_aa[i,j] = sqrt(eps_sq_aa);

                # Get rate factor on central node
                ATT_aa = ATT[i,j];

                # Calculate effective viscosity on ab-nodes
                visc[i,j] = 0.5*(eps_sq_aa)^(p1) * ATT_aa^(p2)

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
function calc_visc_eff_2D_nodes(ux,uy,ATT,H_ice,f_ice,dx,dy;n_glen=3,eps_0=1e-6)

    visc_min = 1e5;

    nx, ny = size(ux);

    # Calculate exponents 
    p1 = (1.0 - n_glen)/(2.0*n_glen);
    p2 = -1.0/n_glen;

    # Calculate squared minimum strain rate 
    eps_0_sq = eps_0*eps_0;

    # Populate strain rates over the whole domain on acx- and acy-nodes

    dudx = fill(0.0,nx,ny);
    dvdy = fill(0.0,nx,ny);
    dudy = fill(0.0,nx,ny);
    dvdx = fill(0.0,nx,ny);
    
    for i = 1:nx
        for j = 1:ny  

            # BC: Periodic boundary conditions
            im1 = i-1
            if im1 == 0
                im1 = nx
            end
            ip1 = i+1
            if ip1 == nx+1
                ip1 = 1
            end

            jm1 = j-1
            if jm1 == 0
                jm1 = ny
            end
            jp1 = j+1
            if jp1 == ny+1
                jp1 = 1
            end

            dudx[i,j] = (ux[ip1,j]-ux[im1,j])/(2.0*dx);
            dudy[i,j] = (ux[i,jp1]-ux[i,jm1])/(2.0*dy);
            dvdx[i,j] = (uy[ip1,j]-uy[im1,j])/(2.0*dx);
            dvdy[i,j] = (uy[i,jp1]-uy[i,jm1])/(2.0*dy);

        end
    end

    # Calculate visc_eff on aa-nodes

    visc   = fill(visc_min,nx,ny);
    eps_aa = fill(eps_0_sq,nx,ny);

    wt0 = 1.0/sqrt(3);
    xn  = [wt0,-wt0,-wt0, wt0];
    yn  = [wt0, wt0,-wt0,-wt0];
    wtn = [1.0,1.0,1.0,1.0];

    for i = 1:nx
        for j = 1:ny  

            if f_ice[i,j] == 1.0

                # BC: Periodic boundary conditions
                im1 = i-1
                if im1 == 0
                    im1 = nx
                end
                ip1 = i+1
                if ip1 == nx+1
                    ip1 = 1
                end

                jm1 = j-1
                if jm1 == 0
                    jm1 = ny
                end
                jp1 = j+1
                if jp1 == ny+1
                    jp1 = 1
                end

                # Get strain rate terms on node locations
                dudxn = acx_to_nodes(dudx,i,j,xn,yn);
                dudyn = acx_to_nodes(dudy,i,j,xn,yn);
                
                dvdxn = acy_to_nodes(dvdx,i,j,xn,yn);
                dvdyn = acy_to_nodes(dvdy,i,j,xn,yn);
                
                # Calculate the total effective strain rate from L19, Eq. 21 
                eps_sq_n = dudxn.^2 + dvdyn.^2 .+ dudxn.*dvdyn .+ 0.25.*(dudyn.+dvdxn).^2 .+ eps_0_sq
                eps_aa = sum(sqrt.(eps_sq_n)) / length(eps_sq_n);

                # Get rate factor on central node
                ATT_aa = ATT[i,j];

                # Calculate effective viscosity on ab-nodes
                viscn = 0.5.*(eps_sq_n).^(p1) .* ATT_aa^(p2);

                visc[i,j] = sum(viscn.*wtn)/sum(wtn);

            end

        end 
    end

    #println("eps: ", extrema(eps_aa))

    return visc

end


function calc_matrix_ssa(ux,uy,N,N_ab,taud_acx,taud_acy,β_acx,β_acy,dx)
    # Calculate the diagnostic SSA velocity solution 
    # given ice thickness, viscosity, driving stress and basal friction coefficient

    # Get some constants 
    nx, ny = size(N);

    dy = dx;
    
    dxdx = (dx*dx);
    dydy = (dy*dy);
    dxdy = (dx*dy);


    #
    #
    # Populate SSA stress balance matrix equation Ax = b 
    # [ A_ux   A_vx ]  [ u ]  = [ b_x ]
    # [ A_uy   A_vy ]  [ v ]    [ b_y ]
    #
    #

    n_terms = 9;
    n_u     = 2*nx*ny;
    n_sprs  = n_u*n_terms;

    # Populate dense array for now
    u = fill(0.0,n_u);
    b = fill(0.0,n_u);

    # Define array A component vectors (I, J, V)
    Ai = fill(0,  n_sprs);
    Aj = fill(0,  n_sprs);
    Av = fill(0.0,n_sprs);

    # Equation is being defined for acx-nodes (x-direction equation)

    k = 0 

    for i in 1:nx
        for j in 1:ny 

            # BC: Periodic boundary conditions
            im1 = i-1
            if im1 == 0
                im1 = nx
            end
            ip1 = i+1
            if ip1 == nx+1
                ip1 = 1
            end

            jm1 = j-1
            if jm1 == 0
                jm1 = ny
            end
            jp1 = j+1
            if jp1 == ny+1
                jp1 = 1
            end

            # Set the row in matrix A that the equation is being defined for:
            nr = (i-1)*ny + j

            # -- vx terms -- 
            
            # ux(i+1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(ip1,j,nx,ny);
            Av[k] = ( 4.0/dxdx*N[ip1,j]);

            # ux(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,j,nx,ny);
            Av[k] = (-4.0/dxdx * (N[ip1,j]+N[i,j])
                     -1.0/dydy * (N_ab[i,j]+N_ab[i,jm1])
                    -β_acx[i,j]);

            # ux(i-1,j)  
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,j,nx,ny);
            Av[k] = ( 4.0/dxdx*N[i,j]);

            # ux(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jp1,nx,ny);
            Av[k] = ( 1.0/dydy*N_ab[i,j]);

            # ux(i,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jm1,nx,ny);
            Av[k] = ( 1.0/dydy*N_ab[i,jm1]);
            
            # -- vy terms -- 

            # uy(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,j,nx,ny);
            Av[k] = (-2.0/dxdy*N[i,j]
                     -1.0/dxdy*N_ab[i,j]);
            
            # uy(i+1,j)
            # NOTE: the N term is `+2`, not `-2` as this reference originally had. With `-2`
            # the operator is neither symmetric nor annihilated by a rigid translation
            # `uy = const`; see the "Linear" section of
            # docs/src/examples/ais-momentum/linear-tuned-gpu.jl. The live port in
            # src/legacy/performance/velocities.jl (`loop1_coeffs`) carries the same fix.
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,j,nx,ny);
            Av[k] = ( 2.0/dxdy*N[ip1,j]
                     +1.0/dxdy*N_ab[i,j]);
            
            # uy(i+1,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,jm1,nx,ny);
            Av[k] = (-2.0/dxdy*N[ip1,j]
                     -1.0/dxdy*N_ab[i,jm1]);
            
            # uy(i,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jm1,nx,ny);
            Av[k] = ( 2.0/dxdy*N[i,j]
                     +1.0/dxdy*N_ab[i,jm1]);

            # [u] value
            u[nr] = ux[i,j];

            # [b] value 
            b[nr] = taud_acx[i,j];

        end
    end

    # Equation is being defined for acy-nodes (y-direction equation)

    for i in 1:nx 
        for j in 1:ny

            # BC: Periodic boundary conditions
            im1 = i-1
            if im1 == 0
                im1 = nx
            end
            ip1 = i+1
            if ip1 == nx+1
                ip1 = 1
            end

            jm1 = j-1
            if jm1 == 0
                jm1 = ny
            end
            jp1 = j+1
            if jp1 == ny+1
                jp1 = 1
            end

            # Set the row in matrix A that the equation is being defined for:
            nr = (i-1)*ny + j + nx*ny

            # -- uy terms -- 
            
            # uy(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jp1,nx,ny);
            Av[k] = ( 4.0/dydy*N[i,jp1]);  
            
            # uy(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,j,nx,ny);
            Av[k] = (-4.0/dydy * (N[i,jp1]+N[i,j])
                     -1.0/dxdx * (N_ab[i,j]+N_ab[im1,j])
                     -β_acy[i,j]);

            # uy(i,j-1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(i,jm1,nx,ny);
            Av[k] = ( 4.0/dydy*N[i,j]);    
            
            # uy(i+1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(ip1,j,nx,ny);
            Av[k] = ( 1.0/dxdx*N_ab[i,j]);     
            
            # uy(i-1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_uy(im1,j,nx,ny);
            Av[k] = ( 1.0/dxdx*N_ab[im1,j]);   
            
            # -- ux terms -- 

            # ux(i,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,jp1,nx,ny);
            Av[k] = ( 2.0/dxdy*N[i,jp1]
                     +1.0/dxdy*N_ab[i,j]);

            # ux(i,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(i,j,nx,ny);
            Av[k] = (-2.0/dxdy*N[i,j]
                     -1.0/dxdy*N_ab[i,j]); 

            # ux(i-1,j+1)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,jp1,nx,ny);
            Av[k] = (-2.0/dxdy*N[i,jp1]
                     -1.0/dxdy*N_ab[im1,j]);
            
            # ux(i-1,j)
            k = k+1;
            Ai[k] = nr;
            Aj[k] = ij2n_ux(im1,j,nx,ny);
            Av[k] = ( 2.0/dxdy*N[i,j]
                     +1.0/dxdy*N_ab[im1,j]);

            # [u] value
            u[nr] = uy[i,j];

            # [b] value 
            b[nr] = taud_acy[i,j];
            
        end
    end
    
    # Now u, b and A components (I, J, V vectors) have been defined.
    # Convert into a sparse array for solving:

    Asp = sparse(Ai,Aj,Av);


    # use_linsolve = true

    # if use_linsolve
    #     prob = LinearProblem(Asp, b; u0=u);
    #     sol = solve(prob);
    #     unew = sol.u;

    # else

    #     unew = Asp \ b;

    # end

    # # Define output velocity arrays with new solution
    # ux1 = fill(0.0,nx,ny);
    # uy1 = fill(0.0,nx,ny);

    # for i = 1:nx
    #     for j = 1:ny
    #         n = ij2n_ux(i,j,nx,ny);
    #         ux1[i,j] = unew[n];
    #         n = ij2n_uy(i,j,nx,ny);
    #         uy1[i,j] = unew[n];
    #     end
    # end

    return Asp, u, b #ux1, uy1
end

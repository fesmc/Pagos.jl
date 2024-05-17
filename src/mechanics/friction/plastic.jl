
"""
    calc_beta_aa_power_plastic(ux_b,uy_b,c_bed,f_ice,q,u_0)

    Calculate basal friction coefficient (beta) that
    enters the SSA solver as a function of basal velocity
    using a power-law form following Bueler and van Pelt (2015)
    
"""
function calc_beta_aa_power_plastic(
    ux_b::Matrix{T},
    uy_b::Matrix{T},
    c_bed::Matrix{T},
    f_ice::Matrix{T},
    q,
    u_0,
) where {T}

    # Local variables
    ub_min = 1e-3               # [m/yr] Small min. velocity > 0 to avoid divide by 0
    ub_sq_min = ub_min^2
    nx, ny = size(ux_b)
    beta = fill(0.0, nx, ny)    # Initially set friction to zero everywhere

    for i = 1:nx
        for j = 1:ny
            im1, jm1 = periodic_bc_minusindices(i, nx), periodic_bc_minusindices(j, ny)

            if f_ice[i, j] == 1.0
                # Fully ice-covered point with some fully ice-covered neighbors 
                cb_aa = c_bed[i, j]

                if q == 1.0
                    # Linear law, no f(ub) term
                    beta[i, j] = cb_aa / u_0
                else
                    # Non-linear law with f(ub) term 
                    # Unstagger velocity components to aa-nodes 
                    ux_aa = 0.5 * (ux_b[i, j] + ux_b[im1, j])
                    uy_aa = 0.5 * (uy_b[i, j] + uy_b[i, jm1])
                    uxy_aa = sqrt(ux_aa^2 + uy_aa^2 + ub_sq_min)

                    if q == 0
                        # Plastic law
                        beta[i, j] = cb_aa * (1.0 / uxy_aa)
                    else
                        beta[i, j] = cb_aa * (uxy_aa / u_0)^q * (1.0 / uxy_aa)
                    end
                end

            else
                # Assign minimum velocity value, no staggering for simplicity

                if q == 1.0
                    # Linear law, no f(ub) term
                    beta[i, j] = c_bed[i, j] / u_0

                else
                    uxy_b = ub_min

                    if q == 0.0
                        # Plastic law
                        beta[i, j] = c_bed[i, j] * (1.0 / uxy_b)
                    else
                        beta[i, j] = c_bed[i, j] * (uxy_b / u_0)^q * (1.0 / uxy_b)
                    end
                end
            end
        end
    end

    return beta
end

function calc_beta_aa_power_plastic_nodes(
    ux_b::Matrix{T},
    uy_b::Matrix{T},
    c_bed::Matrix{T},
    f_ice::Matrix{T},
    q,
    u_0,
) where {T}


    # Local variables
    ub_min = 1e-3               # [m/yr] Minimum velocity is positive small value to avoid divide by zero
    ub_sq_min = ub_min^2

    nx, ny = size(ux_b)

    # Initially set friction to zero everywhere
    beta = fill(0.0, nx, ny)

    wt0 = 1.0 / sqrt(3)
    xn = [wt0, -wt0, -wt0, wt0]
    yn = [wt0, wt0, -wt0, -wt0]
    wtn = [1.0, 1.0, 1.0, 1.0]

    for i = 1:nx
        for j = 1:ny

            if f_ice[i, j] == 1.0
                # Fully ice-covered point with some fully ice-covered neighbors 
                cb_aa = c_bed[i, j]

                if q == 1.0
                    # Linear law, no f(ub) term
                    beta[i, j] = cb_aa / u_0

                else
                    # Non-linear law with f(ub) term 
                    uxn = acx_to_nodes(ux_b, i, j, xn, yn)
                    uyn = acy_to_nodes(uy_b, i, j, xn, yn)
                    uxyn = sqrt.(uxn .^ 2 .+ uyn .^ 2 .+ ub_sq_min)

                    if q == 0
                        # Plastic law
                        betan = cb_aa .* (1.0 ./ uxyn)
                    else
                        betan = cb_aa .* (uxyn ./ u_0) .^ q .* (1.0 ./ uxyn)
                    end
                    beta[i, j] = sum(betan .* wtn) / sum(wtn)
                end

            else
                # Assign minimum velocity value, no staggering for simplicity
                if q == 1.0
                    # Linear law, no f(ub) term
                    beta[i, j] = c_bed[i, j] / u_0
                else
                    uxy_b = ub_min
                    if q == 0.0
                        # Plastic law
                        beta[i, j] = c_bed[i, j] * (1.0 / uxy_b)
                    else
                        beta[i, j] = c_bed[i, j] * (uxy_b / u_0)^q * (1.0 / uxy_b)
                    end
                end
            end
        end
    end

    return beta
end
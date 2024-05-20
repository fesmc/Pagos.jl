
# # Functions to calculate velocity 
# function F_integral!(F, mu, H, f, zeta_aa, m, nx, ny)
#     for i = 1:nx, j = 1:ny
#         F[i, j] = 1 / mu[i, j] * H[i, j]^f * zeta_aa[i, j, k]
#     end
#     return Fn
# end

# function basalvelocity!(state, domain, params, tools)
#     (; H, z_b, ux, uy, β_acx, β_acy) = state
#     (; nx, ny, dx, dy) = domain
#     (; ρ, g) = params
#     (; Ai, Aj, Av, u, b) = tools
#     basalvelocity!(ux, uy, H, z_b, β_acx, β_acy, ρ, g, nx, ny, dx, dy, Ai, Aj, Av, u, b)
#     return nothing
# end

"""
    calc_vel_ssa(ux, uy, H, μ, domain, state)

Calculate the diagnostic SSA velocity solution
given ice thickness, viscosity, driving stress and basal friction coefficient.
"""
function calc_vel_ssa!(icesheet)
    (; state, domain, params, tools) = icesheet
    calc_vel_ssa!(state, domain, params, tools)
    return nothing
end

function calc_vel_ssa!(state, domain, params, tools)
    (; H, μ, taud_acx, taud_acy, β_acx, β_acy, ux, uy) = state
    (; nx, ny, dx, dy) = domain
    (; ρ, g) = params
    (; Ai, Aj, Av, u, b) = tools
    calc_vel_ssa!(ux, uy, H, μ, taud_acx, taud_acy, β_acx, β_acy, dx, dy, Ai, Aj, Av, u, b)
    return nothing
end

function calc_vel_ssa!(ux, uy, H, μ, taud_acx, taud_acy,
    β_acx, β_acy, dx, dy, Ai, Aj, Av, u, b)

    dxdx = (dx * dx)
    dydy = (dy * dy)
    dxdy = (dx * dy)

    # Define vertically-integrated viscosity (aa-nodes)
    N = H .* μ

    # Stagger N to ab-nodes
    N_ab = copy(N)

    @inbounds for i = 1:nx, j = 1:ny
            ip1 = periodic_bc_plusindex(i, nx)
            jp1 = periodic_bc_plusindex(j, ny)
            N_ab[i, j] = 0.25 * (N[i, j] + N[ip1, j] + N[i, jp1] + N[ip1, jp1])
    end

    # Populate SSA stress balance matrix equation Ax = b 
    # [ A_ux   A_vx ]  [ u ]  = [ b_x ]
    # [ A_uy   A_vy ]  [ v ]    [ b_y ]

    # Equation is being defined for acx-nodes (x-direction equation)
    k = 0
    @inbounds for i = 1:nx, j = 1:ny
        im1, ip1, jm1, jp1 = periodic_bc_indices(i, j, nx, ny)

        # Set the row in matrix A that the equation is being defined for:
        nr = (i - 1) * ny + j

        # -- vx terms --

        # ux(i+1,j)
        k += 1
        Ai[k] = nr
        Aj[k] = ij2n_ux(ip1, j, nx, ny)
        Av[k] = (4.0 / dxdx * N[ip1, j])

        # ux(i,j)
        k += 1
        Ai[k] = nr
        Aj[k] = ij2n_ux(i, j, nx, ny)
        Av[k] = (
            -4.0 / dxdx * (N[ip1, j] + N[i, j]) -
            1.0 / dydy * (N_ab[i, j] + N_ab[i, jm1]) - β_acx[i, j]
        )

        # ux(i-1,j)
        k += 1
        Ai[k] = nr
        Aj[k] = ij2n_ux(im1, j, nx, ny)
        Av[k] = (4.0 / dxdx * N[i, j])

        # ux(i,j+1)
        k += 1
        Ai[k] = nr
        Aj[k] = ij2n_ux(i, jp1, nx, ny)
        Av[k] = (1.0 / dydy * N_ab[i, j])

        # ux(i,j-1)
        k += 1
        Ai[k] = nr
        Aj[k] = ij2n_ux(i, jm1, nx, ny)
        Av[k] = (1.0 / dydy * N_ab[i, jm1])

        # -- vy terms --

        # uy(i,j)
        k += 1
        Ai[k] = nr
        Aj[k] = ij2n_uy(i, j, nx, ny)
        Av[k] = (-2.0 / dxdy * N[i, j] - 1.0 / dxdy * N_ab[i, j])

        # uy(i+1,j)
        k += 1
        Ai[k] = nr
        Aj[k] = ij2n_uy(ip1, j, nx, ny)
        Av[k] = (-2.0 / dxdy * N[ip1, j] + 1.0 / dxdy * N_ab[i, j])

        # uy(i+1,j-1)
        k += 1
        Ai[k] = nr
        Aj[k] = ij2n_uy(ip1, jm1, nx, ny)
        Av[k] = (-2.0 / dxdy * N[ip1, j] - 1.0 / dxdy * N_ab[i, jm1])

        # uy(i,j-1)
        k += 1
        Ai[k] = nr
        Aj[k] = ij2n_uy(i, jm1, nx, ny)
        Av[k] = (2.0 / dxdy * N[i, j] + 1.0 / dxdy * N_ab[i, jm1])

        # [u] value
        u[nr] = ux[i, j]

        # [b] value 
        b[nr] = taud_acx[i, j]
    end

    # Equation is being defined for acy-nodes (y-direction equation)
    @inbounds for i = 1:nx, j = 1:ny
            im1, ip1, jm1, jp1 = periodic_bc_indices(i, j, nx, ny)

            # Set the row in matrix A that the equation is being defined for:
            nr = (i - 1) * ny + j + nx * ny

            # -- uy terms -- 

            # uy(i,j+1)
            k += 1
            Ai[k] = nr
            Aj[k] = ij2n_uy(i, jp1, nx, ny)
            Av[k] = (4.0 / dydy * N[i, jp1])

            # uy(i,j)
            k += 1
            Ai[k] = nr
            Aj[k] = ij2n_uy(i, j, nx, ny)
            Av[k] = (
                -4.0 / dydy * (N[i, jp1] + N[i, j]) -
                1.0 / dxdx * (N_ab[i, j] + N_ab[im1, j]) - β_acy[i, j]
            )

            # uy(i,j-1)
            k += 1
            Ai[k] = nr
            Aj[k] = ij2n_uy(i, jm1, nx, ny)
            Av[k] = (4.0 / dydy * N[i, j])

            # uy(i+1,j)
            k += 1
            Ai[k] = nr
            Aj[k] = ij2n_uy(ip1, j, nx, ny)
            Av[k] = (1.0 / dxdx * N_ab[i, j])

            # uy(i-1,j)
            k += 1
            Ai[k] = nr
            Aj[k] = ij2n_uy(im1, j, nx, ny)
            Av[k] = (1.0 / dxdx * N_ab[im1, j])

            # -- ux terms -- 

            # ux(i,j+1)
            k += 1
            Ai[k] = nr
            Aj[k] = ij2n_ux(i, jp1, nx, ny)
            Av[k] = (2.0 / dxdy * N[i, jp1] + 1.0 / dxdy * N_ab[i, j])

            # ux(i,j)
            k += 1
            Ai[k] = nr
            Aj[k] = ij2n_ux(i, j, nx, ny)
            Av[k] = (-2.0 / dxdy * N[i, j] - 1.0 / dxdy * N_ab[i, j])

            # ux(i-1,j+1)
            k += 1
            Ai[k] = nr
            Aj[k] = ij2n_ux(im1, jp1, nx, ny)
            Av[k] = (-2.0 / dxdy * N[i, jp1] - 1.0 / dxdy * N_ab[im1, j])

            # ux(i-1,j)
            k += 1
            Ai[k] = nr
            Aj[k] = ij2n_ux(im1, j, nx, ny)
            Av[k] = (2.0 / dxdy * N[i, j] + 1.0 / dxdy * N_ab[im1, j])

            # [u] value
            u[nr] = uy[i, j]

            # [b] value 
            b[nr] = taud_acy[i, j]
    end

    # Now u, b and A components (I, J, V vectors) have been defined.
    # Convert into a sparse array for solving:
    Asp = sparse(Ai, Aj, Av)
    use_linsolve = true

    if use_linsolve
        prob = LinearProblem(Asp, b; u0 = u)
        sol = solve(prob)
        unew = sol.u
    else
        unew = Asp \ b
    end

    # Define output velocity arrays with new solution
    ux1 = fill(0.0, nx, ny)
    uy1 = fill(0.0, nx, ny)

    for i = 1:nx
        for j = 1:ny
            n = ij2n_ux(i, j, nx, ny)
            ux1[i, j] = unew[n]
            n = ij2n_uy(i, j, nx, ny)
            uy1[i, j] = unew[n]
        end
    end

    return ux1, uy1
end

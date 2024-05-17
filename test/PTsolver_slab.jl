using Pagos

f_eta(beta, H, mu) = beta * H / mu
f_F2(H, mu) = H / (3.0 * mu)
f_ub(rho, g, H, alpha, beta) = rho * g * H * alpha / beta
f_u(rho, g, H, alpha, beta, mu) = ub(rho, g, H, alpha, beta) + rho * g * H^2 * alpha / (3.0 * mu)

struct Slab
    H
    mu
    beta
    alpha
    rho
    g
    eta
    F2
    ub
    u
end

function Slab(; H = 1000.0, mu = 1e5, beta = 1e3, alpha = 1e-3, rho = 910.0, g = 9.81)
    eta = f_eta(beta, H, mu)
    F2 = f_F2(H, mu)
    ub = f_ub(rho, g, H, alpha, beta)
    u = f_u(rho, g, H, alpha, beta, mu)
    return Slab(H, mu, beta, alpha, rho, g, eta, F2, ub, u)
end

function solve_slab(slab, dx; nx = 11, ny = 3)

    lx = nx * dx
    ly = ny * dx
    domain = Domain(Float64, lx, ly, nx, ny)
    state = State(domain)
    params = Params(Float64)
    options = Options(Float64, maxiter = 1000)
    icesheet = IceSheet(state, domain, params, options)
    
    (; state, domain, params, options) = icesheet
    state.H .= slab.H0
    state.mu .= slab.mu0
    state.beta .= slab.beta0
    
    for j = 1:ny
        state.z_b[:, j] = [10000.0 - slab["α"] * (x) for x in domain.x]
    end
    z_srf = state.z_b + state.H
    state.beta_acx, state.beta_acy = stagger_beta(state.beta)

    # Solve for new velocity solution
    ux1, uy1 = calc_vel_ssa(ux, uy, H, μ, taud_acx, taud_acy, β_acx, β_acy, dx)

    println("ux, uy: ", extrema(ux), " | ", extrema(uy))

    return ux1, uy1
end

function plot_var2D(var)

    fig, ax, hm = heatmap(var)
    Colorbar(fig[1, end+1], hm)
    save("test.pdf", fig)

    println("extrema: ", round.(extrema(var), sigdigits = 3))
end

######################################

# Case 1 #
slab1 = Slab()
sol1 = solve_slab(slab1, 5e3)
#plot_var2D(strm1["ux"])

# Case 2 #
slab2 = Slab(H = 500.0, mu = 4e5, beta = 30.0, alpha = 1e-3)
sol2 = solve_slab(slab2, 5e3)
#plot_var2D(strm2["ux"])

######################################
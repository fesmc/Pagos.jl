using Pagos

f_eta(beta, H, mu) = beta * H / mu
f_F2(H, mu) = H / (3.0 * mu)
f_ub(rho, g, H, alpha, beta) = rho * g * H * alpha / beta
f_u(rho, g, H, alpha, beta, mu) = f_ub(rho, g, H, alpha, beta) +
    rho * g * H^2 * alpha / (3.0 * mu)

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

function solve_slab(slab, dx; nx = 11, ny = 3, dtau_scaling = 1.0)

    T = Float64
    lx = nx * dx
    ly = ny * dx
    dy = dx
    domain = Domain(T, lx, ly, dx, dy)
    state = State(domain)
    params = Params{T}()
    options = Options{T}(maxiter = 100_000, printout_every = 1000, dtau_scaling = dtau_scaling)
    icesheet = IceSheet(state, domain, params, options)
    
    (; state, domain, params, options) = icesheet
    state.H .= slab.H
    state.mu .= slab.mu
    state.β .= slab.beta
    for j = 1:ny
        state.z_b[:, j] = [(10000.0 - slab.alpha * x) for x in domain.x]
    end
    pseudo_transient!(icesheet)
    return icesheet
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
sol1 = solve_slab(slab1, 5e3, dtau_scaling = 1.0)
# plot_var2D(strm1["ux"])

# Case 2 #
slab2 = Slab(H = 500.0, mu = 4e5, beta = 30.0, alpha = 1e-3)
sol2 = solve_slab(slab2, 5e3, dtau_scaling = 1e1)
# plot_var2D(strm2["ux"])

######################################
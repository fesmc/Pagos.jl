f_eta(beta, H, mu) = beta * H / mu
f_F2(H, mu) = H / (3.0 * mu)
f_ub(rho, g, H, alpha, beta) = rho * g * H * alpha / beta
f_u(rho, g, H, alpha, beta, mu) = f_ub(rho, g, H, alpha, beta) +
    rho * g * H^2 * alpha / (3.0 * mu)
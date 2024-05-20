
# # Functions to calculate velocity 
# function F_integral!(F, mu, H, f, zeta_aa, m, nx, ny)
#     for i = 1:nx, j = 1:ny
#         F[i, j] = 1 / mu[i, j] * H[i, j]^f * zeta_aa[i, j, k]
#     end
#     return Fn
# end

# function basalvelocity!(state, domain, params, tools)
#     (; H, z_b, ux, uy, beta_acx, beta_acy) = state
#     (; nx, ny, dx, dy) = domain
#     (; ρ, g) = params
#     (; Ai, Aj, Av, u, b) = tools
#     basalvelocity!(ux, uy, H, z_b, beta_acx, beta_acy, ρ, g, nx, ny, dx, dy, Ai, Aj, Av, u, b)
#     return nothing
# end

"""
    stagger_beta!(icesheet)
    stagger_beta!(state, domain)

Stagger the basal friction coefficient to the ab-nodes.
"""
function stagger_beta!(icesheet)
    (; state, domain) = icesheet
    stagger_beta!(state, domain)
    return nothing
end

function stagger_beta!(state, domain)
    (; beta, beta_acx, beta_acy) = state
    (; nx, ny) = domain
    stagger!(beta_acx, beta_acy, beta, nx, ny)
    return nothing
end

"""
    stagger!(outx, outy, in, nx, ny)

Stagger the input field to the ab-nodes.
"""
function stagger!(outx, outy, in, nx, ny)
    for i = 1:nx, j = 1:ny
        ip1 = periodic_bc_plusindex(i, nx)
        jp1 = periodic_bc_plusindex(j, ny)

        outx[i, j] = 0.5 * (in[i, j] + in[ip1, j])
        outy[i, j] = 0.5 * (in[i, j] + in[i, jp1])
    end
    return nothing
end
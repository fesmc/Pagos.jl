
"""
    stagger_beta!(icesheet)
    stagger_beta!(state, domain)
    stagger_beta!(β_acx, β_acy, β, nx, ny)

Stagger the basal friction coefficient to the ab-nodes.
"""
function stagger_beta!(icesheet)
    (; state, domain) = icesheet
    stagger_beta!(state, domain)
    return nothing
end

function stagger_beta!(state, domain)
    (; β, β_acx, β_acy) = state
    (; nx, ny) = domain
    stagger_beta!(β_acx, β_acy, β, nx, ny)
    return nothing
end

function stagger_beta!(β_acx, β_acy, β, nx, ny)
    for i = 1:nx
        for j = 1:ny
            ip1 = periodic_bc_plusindices(i, nx)
            jp1 = periodic_bc_plusindices(j, ny)

            β_acx[i, j] = 0.5 * (β[i, j] + β[ip1, j])
            β_acy[i, j] = 0.5 * (β[i, j] + β[i, jp1])
        end
    end
    return nothing
end
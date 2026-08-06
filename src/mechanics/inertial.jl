oul# @dev TODO: predates the current `IceSheet` layout (references `state`/`domain`/`params`
# fields and a `method` variable that no longer exist); needs to be rewritten against the
# current API before use.
function velocity!(icesheet::IceSheet)
    error("velocity!(::IceSheet) in inertial.jl is not yet implemented for the current IceSheet API")
end

"""
$(TYPEDSIGNATURES)

Compute the time derivative `dudt` of the depth-averaged velocity `u` for an inertial
momentum balance, i.e. the right-hand side of the acceleration equation
`ρ ∂u/∂t = driving stress − basal drag − membrane divergence`. The `IceSheet` method unpacks
the required fields and forwards them to the low-level method. Mutates `dudt` in place.
"""
function inertial_velocity!(dudt, u, t, icesheet::IceSheet)
    # @dev TODO: predates the current `IceSheet` layout (references `domain`/`now`/`params`
    # fields that no longer exist); needs to be rewritten against the current API.
    error("inertial_velocity!(::IceSheet) is not yet implemented for the current IceSheet API")
end

function inertial_velocity!(dudt, u, rho_ice, mu, H, z_s, g, tau_b1, tau_b2, dx1, dx2, mask)

    du1_dt = view(dudt, :, :, 1)
    du2_dt = view(dudt, :, :, 2)
    u1 = view(u, :, :, 1)
    u2 = view(u, :, :, 2)
    du1_dx1 = delx1(u1, dx1, mask)
    du1_dx2 = delx2(u1, dx2, mask)
    du2_dx1 = delx1(u2, dx1, mask)
    du2_dx2 = delx2(u2, dx2, mask)
    rho_ice_inv = 1 / rho_ice
    sec_per_year = 3600 * 24 * 365.25
    c = 0.1

    term11 = rho_ice_inv .* delx1(2 .* mu .* H .* (2 .* du1_dx1 .+ du2_dx2), dx1, mask)
    term12 = rho_ice_inv .* delx2(mu .* H .* (du1_dx2 .+ du2_dx1), dx2, mask)
    term13 = - rho_ice_inv .* tau_b1
    term14 = - g .* H .* delx1(z_s, dx1, mask)
    term15 = - u1 .* H .* delx1(u1, dx1, mask)
    term16 = - u2 .* H .* delx2(u1, dx2, mask)
    println("Extrema: term11 = $(extrema(term11[mask])), term12 = $(extrema(term12[mask])), term13 = $(extrema(term13[mask])), term14 = $(extrema(term14[mask])), term15 = $(extrema(term15[mask])), term16 = $(extrema(term16[mask]))")

    term21 = rho_ice_inv .* delx2(2 .* mu .* H .* (2 .* du2_dx2 .+ du1_dx1), dx2, mask)
    term22 = rho_ice_inv .* delx1(mu .* H .* (du1_dx2 .+ du2_dx1), dx1, mask)
    term23 = - rho_ice_inv .* tau_b2
    term24 = - g .* H .* delx2(z_s, dx2, mask)
    term25 = - u1 .* H .* delx1(u2, dx1, mask)
    term26 = - u2 .* H .* delx2(u2, dx2, mask)
    # @show extrema(term21), extrema(term22), extrema(term23), extrema(term24), extrema(term25), extrema(term26)

    for I in CartesianIndices(du1_dt)
        if H[I] > 0
            du1_dt[I] = (
                term11[I] .* 1e3 +
                term12[I] .* 1e3 +
                term13[I] .* 1e2 +
                term14[I] .* 1e2 +
                term15[I] +
                term16[I]
            ) / H[I]
            du2_dt[I] = (
                term21[I] .* 1e3 +
                term22[I] .* 1e3 +
                term23[I] .* 1e2 +
                term24[I] .* 1e2 +
                term25[I] +
                term26[I]
            ) / H[I]
        else
            du1_dt[I] = 0
            du2_dt[I] = 0
        end
    end

    return nothing

end

function inertial_Hvelocity!(duHdt, u, rho_ice, mu, H, z_s, g, tau_b1, tau_b2, dx1, dx2)

    du1H_dt = view(duHdt, :, :, 1)
    du2H_dt = view(duHdt, :, :, 2)
    u1 = view(u, :, :, 1)
    u2 = view(u, :, :, 2)
    du1_dx1 = delx1(u1, dx1)
    du1_dx2 = delx2(u1, dx2)
    du2_dx1 = delx1(u2, dx1)
    du2_dx2 = delx2(u2, dx2)
    rho_ice_inv = 1 / rho_ice
    sec_per_year = 3600 * 24 * 365.25

    term11 = rho_ice_inv .* delx1(2 .* mu .* H .* (2 .* du1_dx1 .+ du2_dx2), dx1)
    term12 = rho_ice_inv .* delx2(mu .* H .* (du1_dx2 .+ du2_dx1), dx2)
    term13 = - rho_ice_inv .* tau_b1
    term14 = - g .* H .* delx1(z_s, dx1)
    term15 = - u1 .* delx1(u1 .* H, dx1)
    term16 = - u2 .* delx2(u1 .* H, dx2)
    @show extrema(term11), extrema(term12), extrema(term13), extrema(term14), extrema(term15), extrema(term16)
    du1H_dt .= term11 .+
        term12 .+
        term13 .+
        term14 .+
        term15 .+
        term16
    
    term21 = rho_ice_inv .* delx2(2 .* mu .* H .* (2 .* du2_dx2 .+ du1_dx1), dx2)
    term22 = rho_ice_inv .* delx1(mu .* H .* (du1_dx2 .+ du2_dx1), dx1)
    term23 = - rho_ice_inv .* tau_b2
    term24 = - g .* H .* delx2(z_s, dx2)
    term25 = - u1 .* delx1(u2 .* H, dx1)
    term26 = - u2 .* delx2(u2 .* H, dx2)
    @show extrema(term21), extrema(term22), extrema(term23), extrema(term24), extrema(term25), extrema(term26)
    du2H_dt .= term21 .+
        term22 .+
        term23 .+
        term24 .+
        term25 .+
        term26

    # du1H_dt .*= sec_per_year^2
    # du2H_dt .*= sec_per_year^2

    return nothing

end
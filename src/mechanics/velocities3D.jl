"""
$(TYPEDSIGNATURES)

Compute the basal velocity `vb` from the depth-averaged velocity `v`.
"""
function basal_velocity_from_depthavg_velocity!(vb, v, beta, F2)
    @. vb = v / (1 + beta * F2)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Compute the basal velocity `vb` from the surface velocity `vs`.
"""
function basal_velocity_from_surface_velocity!(vb, vs, beta, F1)
    @. vb = vs / (1 + beta * F1)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Compute the surface velocity `vs` from the basal velocity `vb`.
"""
function surface_velocity!(vs, vb, beta, F1)
    @. vs = vb * (1 + beta * F1)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Compute the depth-averaged velocity `v` from the basal velocity `vb`.
"""
function depthavg_velocity!(v, vb, beta, F2)
    @. v = vb * (1 + beta * F2)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Compute the 3D velocity field. `F1` is used as a work array for the running viscosity
integral; it is zeroed on entry. `sigma` must be a host-side vector of sigma-coordinate
midpoints (accessible by index on the CPU during the layer loop).
"""
function velocities3D!(F1, vx3D, vy3D, vb_x, vb_y, mu, beta, H, sigma)
    F1 .= zero(eltype(F1))
    for l in eachindex(sigma)
        aggregate_viscosity_integral!(F1, mu, H, 1, sigma, l)
        layer_velocity!(vx3D, l, vb_x, beta, F1)
        layer_velocity!(vy3D, l, vb_y, beta, F1)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Aggregate the viscosity integral over all sigma layers into `Fm`. `sigma` must be a
host-side vector.
"""
function aggregated_viscosity_integral!(Fm, mu, H, m, sigma)
    Fm .= zero(eltype(Fm))
    for l in eachindex(sigma)
        aggregate_viscosity_integral!(Fm, mu, H, m, sigma, l)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Accumulate the `l`-th layer's contribution to the viscosity integral `Fm` using a
midpoint Riemann sum in sigma coordinates. GPU-compatible.
"""
function aggregate_viscosity_integral!(Fm, mu, H, m, sigma, l)
    dsigma = l == 1 ? sigma[l] : sigma[l] - sigma[l - 1]
    s_minus_z_over_H = 1 - sigma[l] + dsigma / 2
    backend = get_backend(Fm)
    kernel! = _aggregate_viscosity_integral!(backend)
    kernel!(Fm, mu, H, m, s_minus_z_over_H, dsigma, l; ndrange = size(Fm))
    KernelAbstractions.synchronize(backend)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Write the 3D velocity at layer `l` into `v` given the basal velocity `vb` and the
viscosity integral `F1` accumulated up to that layer. GPU-compatible.
"""
function layer_velocity!(v, l, vb, beta, F1)
    backend = get_backend(v)
    kernel! = _layer_velocity!(backend)
    kernel!(v, vb, beta, F1, l; ndrange = size(vb))
    KernelAbstractions.synchronize(backend)
    return nothing
end

@kernel function _aggregate_viscosity_integral!(Fm, mu, H, m, s_minus_z_over_H, dsigma, l)
    i, j = @index(Global, NTuple)
    @inbounds Fm[i, j] += (s_minus_z_over_H ^ m * dsigma * H[i, j]) / mu[i, j, l]
end

@kernel function _layer_velocity!(v, vb, beta, F1, l)
    i, j = @index(Global, NTuple)
    @inbounds v[i, j, l] = vb[i, j] * (1 + beta[i, j] * F1[i, j])
end

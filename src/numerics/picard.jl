
function calc_picard_convergence_metric(
    ux,
    uy,
    ux0,
    uy0;
    mask_acx = nothing,
    mask_acy = nothing,
    vel_tol = 1e-1,
    du_reg = 1e-5,
)

    if isnothing(mask_acx)
        # Assume all points are valid
        kkx = findall(abs.(ux) .> vel_tol)
    else
        kkx = findall(abs.(ux) .> vel_tol .&& mask_acx)
    end

    if isnothing(mask_acy)
        # Assume all points are valid
        kky = findall(abs.(uy) .> vel_tol)
    else
        kky = findall(abs.(uy) .> vel_tol .&& mask_acy)
    end

    res1 = sqrt(sum((ux[kkx] .- ux0[kkx]) .^ 2) + sum((uy[kky] .- uy0[kky]) .^ 2))
    res2 = sqrt(sum((ux[kkx] .+ ux0[kkx]) .^ 2) + sum((uy[kky] .+ uy0[kky]) .^ 2))
    resid = 2.0 * res1 / (res2 + du_reg)

    return resid
end
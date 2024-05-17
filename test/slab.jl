using Printf
using CairoMakie

include(pwd() * "/src/pagos-base.jl");

"""
    define_slab(;H0=1000.0,μ0=1e5,β0=1e3,α=1e-2,ρ=910.0,g=9.81,verbose=true)

    Calculate analytical solution of slab benchmark case
"""
function define_slab(;
    H0 = 1000.0,
    μ0 = 1e5,
    β0 = 1e3,
    α = 1e-2,
    ρ = 910.0,
    g = 9.81,
    verbose = true,
)

    # Get constants
    η = β0 * H0 / μ0
    F2 = H0 / (3.0 * μ0)

    # Calculate rate factor
    # ATT = (2.0*visc_eff)^(-1) 

    # Calculate analytical velocity solution
    ub = ρ * g * H0 * α / β0
    u = ub + ρ * g * H0^2 * α / (3.0 * μ0)

    # Define dictionary with values of analytical stream solution
    strm = Dict(
        "H0" => H0,
        "μ0" => μ0,
        "β0" => β0,
        "α" => α,
        "ρ" => ρ,
        "η" => η,
        "g" => g,
        "F2" => F2,
        "ub" => ub,
        "u" => u,
    )

    # Print summary
    if verbose
        nms = ["H0", "μ0", "β0", "α", "ρ", "ub", "u"]
        println("Analytical stream defined:")
        for nm in nms
            @printf("  %4s = %5.3g", nm, strm[nm])
        end
        print("\n")
    end

    return strm
end

function solve_slab(an, dx; nx = 11, ny = 3)

    # Define axes x,y ###

    xmin = 0.0
    ymin = -(ny - 1) / 2 * dx

    xc = [xmin + (i - 1) * dx for i = 1:nx]
    yc = [ymin + (i - 1) * dx for i = 1:ny]

    # Define variables with initial values ###

    # Ice thickness
    H = fill(an["H0"], nx, ny)

    # Bed elevation
    z_bed = fill(0.0, nx, ny)
    for j = 1:ny
        z_bed[:, j] = [10000.0 - an["α"] * (x) for x in xc]
    end

    # Surface elevation
    z_srf = z_bed .+ H

    # Viscosity
    μ = fill(an["μ0"], nx, ny)

    # Basal friction 
    β = fill(an["β0"], nx, ny)
    β_acx, β_acy = stagger_beta(β)

    # Get driving stress
    taud_acx, taud_acy = calc_driving_stress(H, z_srf, dx, dx, an["ρ"], an["g"])

    # Now, solve for velocity:

    # Define ux, uy and store intial guess 
    ux = fill(0.0, nx, ny)
    uy = fill(0.0, nx, ny)

    # Solve for new velocity solution
    ux1, uy1 = calc_vel_ssa(ux, uy, H, μ, taud_acx, taud_acy, β_acx, β_acy, dx)

    # Update velocity to current solution
    ux, uy = ux1, uy1

    println("ux, uy: ", extrema(ux), " | ", extrema(uy))

    # Collect variables for output 

    return Dict(
        "xc" => xc,
        "yc" => yc,
        "H" => H,
        "z_bed" => z_bed,
        "z_srf" => z_srf,
        "μ" => μ,
        "β" => β,
        "β_acx" => β_acx,
        "β_acy" => β_acy,
        "taud_acx" => taud_acx,
        "taud_acy" => taud_acy,
        "ux" => ux,
        "uy" => uy,
    )
end

function plot_var2D(var)

    fig, ax, hm = heatmap(var)
    Colorbar(fig[1, end+1], hm)
    save("test.pdf", fig)

    println("extrema: ", round.(extrema(var), sigdigits = 3))
end


######################################
# General Parameters #

ρ = 910.0
g = 9.81

######################################

## Test slab ##
if true
    # Case 1 #
    an1 = define_slab(H0 = 1000.0, μ0 = 1e5, β0 = 1e3, α = 1e-3, ρ = ρ, g = g)
    strm1 = solve_slab(an1, 5e3)
    #plot_var2D(strm1["ux"])

    # Case 2 #
    an2 = define_slab(H0 = 500.0, μ0 = 4e5, β0 = 30.0, α = 1e-3, ρ = ρ, g = g)
    strm2 = solve_slab(an2, 5e3)
    #plot_var2D(strm2["ux"])
end
######################################

println("Done.")
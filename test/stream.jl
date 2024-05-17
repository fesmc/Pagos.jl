using NCDatasets
using Printf
using CairoMakie
using Interpolations

include(pwd() * "/src/pagos-base.jl");

"""
    calc_vel_stream!(p)

Calculate analytical solution of stream 
"""
function define_stream_slab(;
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

function solve_stream_slab(an, dx; nx = 11, ny = 3)

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
    ux = fill(0.0, nx, ny)
    uy = fill(0.0, nx, ny)
    ux1, uy1 = calc_vel_ssa(ux, uy, H, μ, taud_acx, taud_acy, β_acx, β_acy, dx)
    ux, uy = ux1, uy1

    #println("ux, uy: ", extrema(ux), " | ", extrema(uy) )

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

function calc_schoof2006_vel(y, taud, H0, B, L, W, m)

    # Calculate the analytical solution for u
    ua = -2.0 * taud^3 * L^4 / (B^3 * H0^3)
    ub = (1.0 / 4.0) * ((y / L)^4.0 - (m + 1.0)^(4.0 / m))
    uc =
        (-3.0 / ((m + 1.0) * (m + 4.0))) *
        (abs(y / L)^(m + 4.0) - (m + 1.0)^(1.0 + (4.0 / m)))
    ud =
        (3.0 / ((m + 1.0)^2 * (2.0 * m + 4.0))) *
        (abs(y / L)^(2 * m + 4.0) - (m + 1.0)^(2.0 + (4.0 / m)))
    ue =
        (-1.0 / ((m + 1.0)^3 * (3.0 * m + 4.0))) *
        (abs(y / L)^(3 * m + 4.0) - (m + 1.0)^(3.0 + (4.0 / m)))
    u = ua * (ub + uc + ud + ue)

    # Outside the ice-stream, velocity is zero
    if abs(y) > W || abs(u) < 1e-10
        u = 0.0
    end

    return u
end

function define_stream_schoof2006(
    dx;
    H0 = 1e3,
    rf = 1e-16,
    α = 1e-3,
    ρ = 910.0,
    g = 9.81,
    n_glen = 3,
    W = 25e3,
    m = 1.55,
)

    # Schoof (2006) domain - constant slope slab

    # Intialize domain 

    # Define axes x,y ###

    ymax = W * 2
    ymin = -ymax

    ny = floor(Int, (ymax - ymin) / dx) + 1
    y0 = -(ny - 1) / 2 * dx
    yc = [ymin + (i - 1) * dx for i = 1:ny]

    # Calculate the gravitational driving stress f
    taud = ρ * g * H0 * α

    # Calculate the ice hardness factor B
    B = rf^(-1 / n_glen)

    # Determine constant L (ice-stream width)
    L = W / ((1.0 + m)^(1.0 / m))

    # Calculate the till yield stress across the stream
    # and analytical velocity solution
    tau_c = taud .* abs.(yc ./ L) .^ m
    ux = calc_schoof2006_vel.(yc, taud, H0, B, L, W, m)

    # Assign analytical values (tau_c as a boundary condition, ux as initial condition)
    c_bed = tau_c

    # Deduce analytical β. 
    β = tau_c ./ ux
    kk = findall(.!isfinite.(β))
    β[kk] .= 0.0

    # Determine constant L too, for diagnostic output
    L = W / ((1.0 + m)^(1.0 / m))

    println("SLAB-S06: H0      = ", H0)
    println("SLAB-S06: alpha   = ", α)
    println("SLAB-S06: W       = ", W)
    println("SLAB-S06: L       = ", L)
    println("SLAB-S06: m       = ", m)
    println("SLAB-S06: ρ g     = ", ρ, " ", g)
    println("SLAB-S06: taud    = ", taud)
    println("SLAB-S06: ATT     = ", rf)
    println("SLAB-S06: c_bed   = ", extrema(c_bed))
    println("SLAB-S06: β       = ", extrema(β))
    println("SLAB-S06: tau_c   = ", extrema(tau_c))
    println("SLAB-S06: ux      = ", extrema(ux))

    return Dict(
        "H0" => H0,
        "α" => α,
        "W" => W,
        "L" => L,
        "m" => m,
        "ρ" => ρ,
        "g" => g,
        "n_glen" => n_glen,
        "rf" => rf,
        "yc" => yc,
        "ux" => ux,
        "tau_c" => tau_c,
        "c_bed" => c_bed,
        "β" => β,
    )
end

function solve_stream_schoof2006(
    an,
    dx;
    xmax = 140e3,
    ymax = maximum(schf_an["yc"]),
    beta_q = 0.0,
    beta_u0 = 1.0,
    eps_0 = 1e-6,
    μ0 = nothing,
    n_iter = 10,
)

    # Intialize domain 

    # Define axes x,y ###

    xmin = 0.0
    ymin = -ymax

    nx = floor(Int, xmax / dx) + 1
    ny = floor(Int, (ymax - ymin) / dx) + 1

    x0 = 0.0
    y0 = -(ny - 1) / 2 * dx

    xc = [xmin + (i - 1) * dx for i = 1:nx]
    yc = [ymin + (i - 1) * dx for i = 1:ny]

    # ===== Intialize topography and set parameters =========

    H = fill(an["H0"], nx, ny)

    z_bed = fill(NaN, nx, ny)

    for i = 1:nx
        z_bed[i, :] .= 10000.0 .- an["α"] .* (xc[i] .- xmin)
    end

    # Surface elevation
    z_srf = z_bed .+ H

    # Ice fraction mask 
    f_ice = fill(1.0, nx, ny)

    # Driving stress
    taud_acx, taud_acy = calc_driving_stress(H, z_srf, dx, dx, an["ρ"], an["g"])

    # Bed properties 

    interp_linear = linear_interpolation(an["yc"], an["c_bed"])

    c_bed = fill(0.0, nx, ny)

    i = 1
    c_bed[1, :] = interp_linear.(yc)

    for i = 2:nx
        c_bed[i, :] .= c_bed[1, :]
    end

    # Rate factor 
    ATT = fill(an["rf"], nx, ny)

    # Viscosity
    if isnothing(μ0)
        const_visc = false
    else
        const_visc = true
    end

    # Now, solve for velocity:
    μ = fill(0.0, nx, ny)
    β = fill(0.0, nx, ny)
    β_acx = fill(0.0, nx, ny)
    β_acy = fill(0.0, nx, ny)
    ux = fill(0.0, nx, ny)
    uy = fill(0.0, nx, ny)

    # Initialize solution to analytical result 

    interp_linear = linear_interpolation(an["yc"], an["ux"])

    i = 1
    ux[1, :] = interp_linear.(yc)

    for i = 2:nx
        ux[i, :] .= ux[1, :]
    end

    # Relaxation weighting
    f_rel = 0.70
    f_rel_visc = 0.7

    # Tolerance for stopping iterations 
    iter_tol = 1e-3

    for iter = 1:n_iter

        # Store solution from previous iteration
        ux0, uy0 = ux, uy

        # Calculate viscosity
        if const_visc
            μ .= μ0
        else
            μ_new = calc_visc_eff_2D_aa(
                ux,
                uy,
                ATT,
                H,
                f_ice,
                dx,
                dx;
                n_glen = an["n_glen"],
                eps_0 = eps_0,
            )
            #μ_new = calc_visc_eff_2D_nodes(ux,uy,ATT,H,f_ice,dx,dx;n_glen=an["n_glen"],eps_0=eps_0);
            #μ = f_rel_visc .* μ_new + (1-f_rel) .* μ;
            μ = μ_new
        end

        # Calculate basal friction
        β_new = calc_beta_aa_power_plastic(ux, uy, c_bed, f_ice, beta_q, beta_u0)
        #β_new = calc_beta_aa_power_plastic_nodes(ux,uy,c_bed,f_ice,beta_q,beta_u0);
        #β = f_rel .* β_new + (1-f_rel) .* β;
        β = β_new
        β_acx, β_acy = stagger_beta(β)

        # Calculate new velocity solution
        ux_new, uy_new =
            calc_vel_ssa(ux, uy, H, μ, taud_acx, taud_acy, β_acx, β_acy, dx)

        # Relax towards new solution 
        ux = f_rel .* ux_new .+ (1 - f_rel) .* ux
        uy = f_rel .* uy_new .+ (1 - f_rel) .* uy

        # Check convergence
        #du = sqrt( sum((ux-ux0).^2 + (uy-uy0).^2) / (nx*ny) );
        du = calc_picard_convergence_metric(ux, uy, ux0, uy0)

        @printf(
            "iter: %5i %8.3e | %8.3e %8.3e | %8.3e %8.3e\n",
            iter,
            du,
            extrema(ux)...,
            extrema(β)...
        )
        #println("iter: ",iter, " ", du, " ", extrema(ux), " | ", extrema(β))

        if du <= iter_tol && iter > 1
            break
        end
    end

    # Get final basal stress field 
    taub_acx = ux .* β_acx
    taub_acy = uy .* β_acy
    taub = fill(0.0, nx, ny)
    for i = 1:nx
        for j = 1:ny

            # BC: Periodic boundary conditions
            im1 = i - 1
            if im1 == 0
                im1 = nx
            end
            ip1 = i + 1
            if ip1 == nx + 1
                ip1 = 1
            end

            jm1 = j - 1
            if jm1 == 0
                jm1 = ny
            end
            jp1 = j + 1
            if jp1 == ny + 1
                jp1 = 1
            end

            taub[i, j] = sqrt(
                0.5 * (taub_acx[i, j] + taub_acx[im1, j])^2 +
                0.5 * (taub_acy[i, j] + taub_acy[i, jm1])^2,
            )
        end
    end

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
        "taub_acx" => taub_acx,
        "taub_acy" => taub_acy,
        "taub" => taub,
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
if false
    # Case 1 #
    an1 = define_stream_slab(H0 = 1000.0, μ0 = 1e5, β0 = 1e3, α = 1e-3, ρ = ρ, g = g)
    strm1 = solve_stream_slab(an1, 5e3)
    plot_var2D(strm1["ux"])

    # Case 2 #
    an2 = define_stream_slab(H0 = 500.0, μ0 = 4e5, β0 = 30.0, α = 1e-3, ρ = ρ, g = g)
    strm2 = solve_stream_slab(an2, 5e3)
    plot_var2D(strm2["ux"])
end
######################################

begin
    ## Test schoof2006 ##

    # Get analytical solution

    # Berends et al (2022) IMAU-ICE parameters
    # The solver converges easily to this set of parameters
    schf_an = define_stream_schoof2006(
        1e3;
        H0 = 2e3,
        rf = 1e-18,
        α = 3e-4,
        ρ = 910.0,
        g = 9.81,
        n_glen = 3,
        W = 300e3,
        m = 1.0,
    )

    # Lipscomb et al (2019) CISM parameters
    # The solver does not converge easily to this set of parameters
    #schf_an = define_stream_schoof2006(1e3;H0=1e3,rf=1e-16,α=1e-3,ρ=910.0,g=9.81,n_glen=3,W=25e3,m=1.55);

    # Get numerical solutions

    dxs = [16e3, 8e3, 4e3]
    schf = []

    for k in eachindex(dxs)
        dx = dxs[k]
        tmp = solve_stream_schoof2006(schf_an, dx; n_iter = 200, eps_0 = 1e-6)
        push!(schf, tmp)
    end


    ### Initialize figure ###

    fig = Figure(resolution = (1000, 600))

    lwds = [4, 3, 2.5]
    cols = [:pink, :red, :darkred]

    ## Panel 1: Velocity solution ## 

    ax1 = Axis(fig[1, 1], xlabel = "y (km)", ylabel = "x-velocity (m/yr)")
    xlims!(ax1, extrema(schf_an["yc"]) .* 1e-3)
    #ylims!(ax1,(0,1000))
    ylims!(ax1, extrema(schf_an["ux"]) .* 1.1)
    #ax1.yticks=0:100:2000;

    # Analytical solution
    lines!(
        ax1,
        schf_an["yc"] * 1e-3,
        schf_an["ux"],
        color = :grey60,
        linewidth = 8,
        label = "Analytical solution",
    )

    # Numerical solutions
    for k in eachindex(dxs)
        dx = dxs[k] * 1e-3
        i = floor(Int, length(schf[k]["xc"]) / 2)
        lines!(
            ax1,
            schf[k]["yc"] * 1e-3,
            schf[k]["ux"][i, :],
            color = cols[k],
            linewidth = lwds[k],
            label = "dx = $dx km",
        )
    end

    ## Panel 2: Basal stress ##

    ax2 = Axis(fig[2, 1], xlabel = "y (km)", ylabel = "Basal stress (kPa)")
    ylims!(ax2, (0, 70))

    # Analytical solution
    lines!(
        ax2,
        schf_an["yc"] * 1e-3,
        schf_an["tau_c"] * 1e-3,
        color = :grey60,
        linewidth = 8,
        label = "Analytical solution",
    )

    # Numerical solutions
    for k in eachindex(dxs)
        dx = dxs[k] * 1e-3
        i = floor(Int, length(schf[k]["xc"]) / 2)
        lines!(
            ax2,
            schf[k]["yc"] * 1e-3,
            schf[k]["taub_acx"][i, :] * 1e-3,
            color = cols[k],
            linewidth = lwds[k],
            label = "dx = $dx km",
        )
    end

    axislegend(ax2; position = :ct)

    ## Panel 3: Viscosity ##
    ax3 =
        Axis(fig[1, 2], xlabel = "y (km)", ylabel = "Viscosity (Pa yr)", yscale = log10)

    # Numerical solutions
    for k in eachindex(dxs)
        dx = dxs[k] * 1e-3
        i = floor(Int, length(schf[k]["xc"]) / 2)
        lines!(
            ax3,
            schf[k]["yc"] * 1e-3,
            schf[k]["μ"][i, :],
            color = cols[k],
            linewidth = lwds[k],
            label = "dx = $dx km",
        )
    end

    save("test.pdf", fig)
end



println("Done.")
module PagosMakieExt

using Pagos, Makie, LinearAlgebra
using DocStringExtensions

function Pagos.plot_rate_factor(T, A)
    set_theme!(theme_latexfonts())
    fig = Figure(size = (500, 400))
    ax = Axis(fig[1, 1],
        xlabel = L"Temperature $T'$ (°C)",
        ylabel = L"Rate Factor A (Pa$^{-3}$ s$^{-1}$)",
        yscale = log10,
        title = "Arrhenius Rate Factor for Ice Viscosity"
    )
    lines!(ax, T, A, linewidth = 2)
    ax.xticks = -50:10:50
    ax.yminorticks = IntervalsBetween(10)
    ax.xminorgridvisible = true
    ax.yminorgridvisible = true
    xlims!(ax, -50, 0)
    ylims!(ax, 1e-27, 1e-23)
    return fig
end

function Pagos.plot_melting_point(p, Tm)
    set_theme!(theme_latexfonts())
    fig = Figure(size = (500, 400))
    ax = Axis(fig[1, 1],
        xlabel = L"Pressure $p$ (MPa)",
        ylabel = L"Melting Point $T_m$ (K)",
        title = "Pressure Melting Point of Ice"
    )
    lines!(ax, p ./ 1f6, Tm, linewidth = 2)
    return fig
end

function Pagos.plot_ice_viscosity(η_ice, σ_e, T;
    rate_factor_str = "Arrhenius",
    flow_law_str = "Regularized Glen-Nye",
)
    set_theme!(theme_latexfonts())
    fig = Figure(size = (500, 400))
    ax = Axis(fig[1, 1];
        xlabel = L"Effective stress $\sigma_\mathrm{e}$ (kPa)",
        ylabel = L"Ice viscosity $\eta$ (Pa s)",
        yscale = log10,
        title = "$rate_factor_str Rate Factor & $flow_law_str flow law"
    )
    for i in eachindex(T)
        lines!(ax, σ_e ./ 1f3, η_ice[i], label = "T = $(T[i]) °C")
    end
    ax.xticks = 0:20:100
    ax.yminorticks = IntervalsBetween(10)
    ax.xminorgridvisible = true
    ax.yminorgridvisible = true
    xlims!(ax, 0, 100)
    ylims!(ax, 1e13, 1e17)
    axislegend(ax, position = :rt)

    return fig
end

function Pagos.plot_basal_shear_stress(v_basal, τ_basal, labels)
    set_theme!(theme_latexfonts())
    fig = Figure(size = (600, 400))
    ax = Axis(fig[1, 1];
        xlabel = L"Basal sliding speed $v_b$ (m/a)",
        ylabel = L"Basal shear stress $\tau_b$ (kPa)",
        title = "Basal shear stress vs basal sliding speed",
    )
    for i in eachindex(τ_basal)
        lines!(ax, [norm(v) for v in v_basal], [norm(τ) for τ in τ_basal[i]], label = labels[i])
    end
    axislegend(ax, position = :rb)
    return fig
end

end     # module PagosMakieExt
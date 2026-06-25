"""
$(TYPEDSIGNATURES)

Compute the effective ice thickness ``H_\mathrm{eff}``:

```math
\begin{aligned}
H_\mathrm{eff} = \begin{cases}\dfrac{H_\mathrm{ice}}{f_\mathrm{ice}} & \text{for } f_\mathrm{ice} > 0 \\ 0 & \text{for } f_\mathrm{ice} = 0 \end{cases}
\end{aligned}
```
"""
function H_ice_effective(H_ice, f_ice)
    if f_ice > 0
        H_eff = H_ice / f_ice
    else
        H_eff = 0
    end
    return H_eff
end

"""
($(TYPEDSIGNATURES)
Compute the height of the ice above floatation:

```math
\begin{aligned}
H_\mathrm{grnd} = \max\left(H_\mathrm{ice} - \dfrac{\rho_\mathrm{seawater}}{\rho_\mathrm{ice}} \\max(z_\mathrm{sl} - z_\mathrm{bed}, 0), 0\right)
\end{aligned}
```
"""
function height_above_floatation(H, ρ_ice, ρ_seawater, z_sl, z_bed)
    return max(H - ρ_seawater / ρ_ice * max(z_sl - z_bed, 0), 0)
end

function surface_elevation(H_grnd, H_eff, ρ_ice, ρ_seawater, z_sl, z_bed)
    if H_grnd > 0
        z_srf = z_bed + H_eff
    else
        z_srf = z_sl + (1 - ρ_seawater / ρ_ice) * H_eff
    end
    return z_srf
end

function maximal_surface_elevation(ρ_ice, ρ_seawater, z_sl, z_bed)
    return max(z_bed + H_eff, z_sl + (1 - ρ_seawater / ρ_ice) * H_eff)
end

function H_ice_grounded(H_eff, z_sl, z_bed, ρ_ice, ρ_seawater)
    H_grnd = height_above_floatation(H_eff, ρ_ice, ρ_seawater, z_sl, z_bed)
    return H_grnd
end
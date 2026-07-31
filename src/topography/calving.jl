##########################################################################
# Structs
##########################################################################

"""
$(TYPEDSIGNATURES)

An abstract type to dispatch calving laws via [`calving_rate`](@ref).
"""
abstract type AbstractCalving end

"""
$(TYPEDSIGNATURES)

Struct to specify a constant calving rate:

```math
\\begin{aligned}
\\dot{c} = \\mathrm{const.}
\\end{aligned}
```

# Fields
 - `rate::T`: Prescribed calving rate
"""
@kwdef struct PrescribedCalving{T} <: AbstractCalving
    rate::T = 0.0           # m yr-1
end

"""
$(TYPEDSIGNATURES)

Calving law based on a relaxation `timescale`:

```math
\\begin{aligned}
\\dot{c} = -\\dfrac{H}{\\tau_{\\mathrm{c}}} \\quad \\text{if } H > H_{\\mathrm{crit}}
\\end{aligned}
```

# Fields
 - `H_critical::T`: Critical thickness ``H_{\\mathrm{crit}}`` for calving
 - `timescale::T`: Calving timescale ``\\tau_{\\mathrm{c}}``
 - `max_rate::T`: Maximum calving rate
"""
@kwdef struct RelaxedCalving{T} <: AbstractCalving
    H_critical::T = 0.0         # m
    timescale::T = 1.0          # yr
    max_rate::T = Inf           # m yr-1
end

"""
$(TYPEDSIGNATURES)

Calving law based on a threshold thickness `H_critical` and a `timescale` as in Eq. 3 of [robinson_description_2020](@citet):

```math
\\begin{aligned}
\\dot{c} = -\\dfrac{H_{\\mathrm{eff}} - H_{\\mathrm{crit}}}{\\tau} \\quad \\text{if } H_{\\mathrm{eff}} > H_{\\mathrm{crit}}
\\end{aligned}
```

# Fields
 - `H_critical::T`: Reference thickness
 - `timescale::T`: Calving time scale
 - `max_rate::T`: Maximum calving rate
"""
@kwdef struct ThicknessCalving{T} <: AbstractCalving
    H_critical::T = 100.0       # m
    timescale::T = 1.0          # yr
    max_rate::T = Inf           # m yr-1
end

"""
$(TYPEDSIGNATURES)

Calving law that removes all ice at or below flotation:
```math
\\begin{aligned}
\\dot{c} = -\\dfrac{H}{\\tau_{\\mathrm{c}}} \\quad \\text{if } H_{\\mathrm{eff}} \\leq 0
\\end{aligned}
```

# Fields
 - `timescale::T`: Calving timescale (``\\mathrm{yr}``).
 - `max_rate::T`: Maximum calving rate (``\\mathrm{m}\\,\\mathrm{yr}^{-1}``).
"""
@kwdef struct FlotationCalving{T} <: AbstractCalving
    timescale::T = 1.0       # yr
    max_rate::T = Inf        # m yr-1
end

"""
$(TYPEDSIGNATURES)

Calving law based on the Von Mises stress criterion as in [lipscomb_description_2019](@citet), Eq. 73-75:

```math
\\begin{aligned}
\\dot{c} = -\\dfrac{H_{\\mathrm{eff}} \\, k_{\\tau} \\, \\tau_{\\mathrm{eff}}}{\\sqrt{\\Delta x \\, \\Delta y}}
\\end{aligned}
```

with ``\\tau_{\\mathrm{eff}}`` the effective calving stress, and ``k_{\\tau}`` an empirical calving coefficient. The former is computed as:

```math
\\begin{aligned}
\\tau_{\\mathrm{eff}} = \\max(\\tau_1, 0)^2 + w_2 \\, \\max(\\tau_2, 0)^2
\\end{aligned}
```

# Fields
 - `k_τ::T`: Empirical calving coefficient
 - `w2::T`: Empirical stress coefficient
 - `timescale::T`: Calving time scale
 - `max_rate::T`: Maximum calving rate
"""
@kwdef struct LipscombCalving{T} <: AbstractCalving
    k_τ::T = 0.0025             # m yr-1 Pa-1
    w2::T = 25                  # 1
    timescale::T = 1.0          # yr
    max_rate::T = Inf
end

"""
$(TYPEDSIGNATURES)

Calving law based on the principal strain rate `ε_eff` as in [levermann_kinematic_2012](@citet), Eq. 2:
```math
\\begin{aligned}
\\dot{c} = -\\dfrac{H_{\\mathrm{eff}} \\, k_2 \\, \\varepsilon_{\\mathrm{eff}}}{\\sqrt{\\Delta x \\, \\Delta y}}
\\end{aligned}
```

# Fields
 - `k2::T`: Empirical calving coefficient
 - `max_rate::T`: Maximum calving rate
"""
@kwdef struct LevermannCalving{T} <: AbstractCalving
    k2::T = 0.00025             # m yr-1 Pa-1
    max_rate::T = Inf           # m yr-1
end

"""
$(TYPEDSIGNATURES)

Calving law based on the ice thickness above sea level `H_c` as in [crawford_marine_2021](@citet), Eq. 1:
```math
\\begin{aligned}
\\dot{c} = -I \\, H_{\\mathrm{c}}^{\\alpha}
\\quad \\text{if } H_{\\mathrm{c}} > H_{\\mathrm{crit}} \\text{ and } z_{\\mathrm{bed}} < z_{\\mathrm{sl}}
\\end{aligned}
```

where ``H_{\\mathrm{c}} = z_{\\mathrm{srf}} - z_{\\mathrm{sl}}`` is the subaerial cliff height.
The large exponent ``\\alpha \\approx 7.3`` encodes the structural fragility of tall
ice cliffs: calving is negligible for modest cliff heights but grows explosively once
``H_{\\mathrm{c}}`` exceeds the critical threshold.

# Fields
 - `I::T`: Calving coefficient (``\\mathrm{m}^{1-\\alpha}\\,\\mathrm{yr}^{-1}``).
 - `α::T`: Power-law exponent (dimensionless).
 - `H_critical::T`: Critical subaerial cliff height ``H_{\\mathrm{crit}}`` (``\\mathrm{m}``).
 - `max_rate::T`: Maximum calving rate (``\\mathrm{m}\\,\\mathrm{yr}^{-1}``).

!!! warning
    This does not depend on the flow regime for now and only represents an upper bound on MICI.
"""
@kwdef struct CrawfordCalving{T} <: AbstractCalving
    I::T = 1.9e-16              # m
    α::T = 7.3                  # 1
    H_critical::T = 135.0         # m
    max_rate::T = Inf           # m yr-1
end

"""
$(TYPEDSIGNATURES)

Calving law based on the formulation by [bassis_upper_2012](@citet):
```math
\\begin{aligned}
\\dot{c} = \\dfrac{f_{\\mathrm{ice}} \\, \\max(H_{\\mathrm{eff}} - H_{\\mathrm{max}}, 0)}{\\tau_{\\mathrm{c}}}
\\end{aligned}
```

# Fields
 - `C0::T`: Cohesive strength parameter
 - `α::T`: Friction parameter
 - `r::T`: Crevasse parameter
 - `max_rate::T`: Maximum calving rate
"""
@kwdef struct BassisCalving{T} <: AbstractCalving
    C0::T = 1f6                 # Pa
    α::T = 0.0                  # 1
    r::T = 0.0                  # 1
    max_rate::T = Inf           # m yr-1
end

"""
$(TYPEDSIGNATURES)

Calving law based on the eigencalving criterion by [winkelmann_analytical_2011](@citet):
```math
\\begin{aligned}
\\dot{c} = K \\, \\max(\\dot{\\varepsilon}_1, 0) \\, \\max(\\dot{\\varepsilon}_2, 0)
\\end{aligned}
```
where ``\\dot{\\varepsilon}_1, \\dot{\\varepsilon}_2`` are the two principal horizontal strain rates.
The calving rate vanishes wherever either principal strain rate is compressive.

# Fields
 - `K::T`: Eigencalving coefficient (``\\mathrm{m}\\,\\mathrm{yr}``).
 - `max_rate::T`: Maximum calving rate (``\\mathrm{m}\\,\\mathrm{yr}^{-1}``).
"""
@kwdef struct EigenCalving{T} <: AbstractCalving
    K::T = 1e7               # m yr
    max_rate::T = Inf        # m yr-1
end


"""
$(TYPEDSIGNATURES)

Calving law following [deconto_contribution_2016](@citet)
and [pollard_potential_2015](@citet). When the subaerial cliff height
``H_{\\mathrm{s}} = z_{\\mathrm{srf}} - z_{\\mathrm{sl}}`` exceeds the critical threshold ``H_{\\mathrm{c}}`` (set by the
structural yield strength of ice, \\approx 90\\text{--}100\\,\\text{m}), the cliff fails and
calves at a rate proportional to the excess height:

```math
\\begin{aligned}
\\dot{c} = -\\dfrac{\\max(H_{\\mathrm{s}} - H_{\\mathrm{c}},\\, 0)}{\\tau_{\\mathrm{c}}}
\\quad \\text{if } z_{\\mathrm{bed}} < z_{\\mathrm{sl}}
\\end{aligned}
```

# Fields
 - `H_c::T`: Critical cliff height above sea level (``\\mathrm{m}``).
 - `timescale::T`: Calving timescale (``\\mathrm{yr}``).
 - `max_rate::T`: Maximum calving rate (``\\mathrm{m}\\,\\mathrm{yr}^{-1}``).
"""
@kwdef struct PollardDeContoCalving{T} <: AbstractCalving
    H_c::T = 100.0           # m
    timescale::T = 1.0       # yr
    max_rate::T = Inf        # m yr-1
end

############################################################################
# Functions
############################################################################

"""
$(TYPEDSIGNATURES)

Calculate the calving flux based on the ice thickness `H` and the calving law `calving<:AbstractCalving`.
"""
function calving_rate(calving::PrescribedCalving)
    return -calving.rate
end

function calving_rate(H, calving::RelaxedCalving)
    (; H_critical, timescale, max_rate) = calving
    if H > H_critical
        rate = H / timescale
    else
        rate = 0.0
    end
    return -saturate(rate, 0, max_rate)
end

function calving_rate(H, H_eff, calving::ThicknessCalving)
    (; H_critical, timescale, max_rate) = calving
    H_eff_dt = (H_eff - H_critical) ./ timescale
    rate = min(H / H_eff, 1) * H_eff_dt
    return -saturate(rate, 0, max_rate)
end

function calving_rate(H, H_eff, tau1, tau2, dx, dy, calving::LipscombCalving)
    (; k_τ, w2, max_rate) = calving
    tau_eff = sqrt(max(tau1, 0)^2 + w2 * max(tau2, 0)^2)
    rate = H_eff * k_τ * tau_eff / sqrt(dx * dy)
    return -saturate(rate, 0, max_rate)
end

function calving_rate(H_eff, eps_eff, dx, dy, calving::LevermannCalving)
    (; k2, max_rate) = calving
    calving_ref = max(k2 * eps_eff, 0)
    rate = H_eff * calving_ref / sqrt(dx * dy)
    return -saturate(rate, 0, max_rate)
end

function calving_rate(z_srf, z_sl, z_bed, calving::CrawfordCalving)
    (; I, α, H_critical, max_rate) = calving
    H_c = z_srf - z_sl
    if H_c > H_critical && z_bed < z_sl
        c_dt = I * H_c ^ α
    else
        c_dt = 0.0
    end
    return -saturate(c_dt, 0, max_rate)
end

# TODO: this needs a `seawater_depth` implementation, and the formula re-checked, before
# it can be used (see roadmaps/todo.md).
function calving_rate(H_eff, z_sl, z_bed, f_ice, c, calving::BassisCalving)
    error("calving_rate(::BassisCalving) is not yet implemented: seawater_depth is undefined")

    # (; C0, α, r, max_rate) = calving
    # (; ρ_ice, ρ_seawater, g) = c

    # # TODO: this needs to be checked and seawater depth should be removed
    # H_ocn = seawater_depth(ρ_seawater / ρ_ice, z_sl - z_bed)
    # timescale = C0 + 0.5 * α * ρ_ice * g * H_eff
    # H_max = (1 - r) * timescale / ρ_ice * g +
    #     sqrt((1 - r)^2 * (timescale / ρ_ice * g)^2 + ρ_seawater / ρ_ice * H_ocn ^ 2)

    # if H_eff <= H_max
    #     return 0.0
    # else
    #     return f_ice * max(H_eff - H_max, 0) / timescale
    # end
end

function calving_rate(eps1, eps2, calving::EigenCalving)
    (; K, max_rate) = calving
    rate = K * max(eps1, 0.0) * max(eps2, 0.0)
    return -saturate(rate, 0.0, max_rate)
end

function calving_rate(H, H_eff, calving::FlotationCalving)
    (; timescale, max_rate) = calving
    rate = H_eff <= 0 ? H / timescale : 0.0
    return -saturate(rate, 0.0, max_rate)
end

# TODO: this should simply take the current state as input instead of the individual fields. This however requires to have a clear definition of the state.
"""
$(TYPEDSIGNATURES)

Update the calving flux `c_dt` based on the calving law `calving<:AbstractCalving` and the ice thickness `H`.
"""
function calving_rate!(c_dt, calving::PrescribedCalving)
    pointwise!(calving_rate, c_dt, (), (calving,))
    return nothing
end

function calving_rate!(c_dt, H, calving::RelaxedCalving)
    pointwise!(calving_rate, c_dt, (H,), (calving,))
    return nothing
end

function calving_rate!(c_dt, H, H_eff, calving::ThicknessCalving)
    pointwise!(calving_rate, c_dt, (H, H_eff), (calving,))
    return nothing
end

function calving_rate!(c_dt, H, H_eff, tau1, tau2, dx, dy, calving::LipscombCalving)
    pointwise!(calving_rate, c_dt, (H, H_eff, tau1, tau2), (dx, dy, calving))
    return nothing
end

function calving_rate!(c_dt, H_eff, eps_eff, dx, dy, calving::LevermannCalving)
    pointwise!(calving_rate, c_dt, (H_eff, eps_eff), (dx, dy, calving))
    return nothing
end

function calving_rate!(c_dt, z_srf, z_sl, z_bed, calving::CrawfordCalving)
    pointwise!(calving_rate, c_dt, (z_srf, z_sl, z_bed), (calving,))
    return nothing
end

function calving_rate!(c_dt, H_eff, z_sl, z_bed, f_ice, c, calving::BassisCalving)
    pointwise!(calving_rate, c_dt, (H_eff, z_sl, z_bed, f_ice), (c, calving))
    return nothing
end

function calving_rate!(c_dt, eps1, eps2, calving::EigenCalving)
    pointwise!(calving_rate, c_dt, (eps1, eps2), (calving,))
    return nothing
end

function calving_rate!(c_dt, H, H_eff, calving::FlotationCalving)
    pointwise!(calving_rate, c_dt, (H, H_eff), (calving,))
    return nothing
end

function calving_rate(z_srf, z_sl, z_bed, calving::PollardDeContoCalving)
    (; H_c, timescale, max_rate) = calving
    H_s = z_srf - z_sl
    if H_s > H_c && z_bed < z_sl
        rate = (H_s - H_c) / timescale
    else
        rate = 0.0
    end
    return -saturate(rate, 0.0, max_rate)
end

function calving_rate!(c_dt, z_srf, z_sl, z_bed, calving::PollardDeContoCalving)
    pointwise!(calving_rate, c_dt, (z_srf, z_sl, z_bed), (calving,))
    return nothing
end
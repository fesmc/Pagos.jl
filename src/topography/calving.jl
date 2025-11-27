############################################################################
# Structs
############################################################################

"""
$(TYPEDSIGNATURES)

An abstract type to dispatch calving laws via [`calving_rate`](@ref).
"""
abstract type AbstractCalving end

# TODO: Decide whether we should use this distinction
abstract type AbstractGroundedCalving <: AbstractCalving end
abstract type AbstractFloatingCalving <: AbstractCalving end
"""
$(TYPEDSIGNATURES)

Struct to specify a constant calving rate:

```math
\\begin{aligned}
\\dot{c} = \\mathrm{const.}
\\end{aligned}
```

# Fields
 - `rate::T`: Constant calving rate
"""
@kwdef struct ConstantCalving{T} <: AbstractCalving
    rate::T = 0.0           # m yr-1
end

"""
$(TYPEDSIGNATURES)

Calving law based on a relaxation `timescale`:

```math
\\begin{aligned}
\\dot{c} = -\\dfrac{H}{\\tau_c}
\\end{aligned}
```

# Fields
 - `H_critical::T`: Critical thickness for calving
 - `timescale::T`: Calving time scale
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
\\dot{c} = -\\dfrac{H - H_{ref}}{\\tau}
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

Calving law based on the Von Mises stress criterion as in [lipscomb_description_2019](@citet), Eq. 73-75:
```math
\\begin{aligned}
\\dot{c} = -\\dfrac{H_{eff} k_{\\tau} \\tau_{eff}}{\\sqrt{\\Delta x \\Delta y}}
\\end{aligned}
```

with `\\tau_{eff}` the effective calving stress, and `k_{\\tau}` an empirical calving coefficient. The former is computed as:

```math
\\begin{aligned}
\\tau_{eff} = max(\\tau_1, 0)^2 + w_2 max(\\tau_2, 0)^2
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

Calving law based on the principal strain rate `ε_eff` as in [levermann_calving_2012](@citet), Eq. 2:
```math
\\begin{aligned}
\\dot{c} = -\\dfrac{H_{eff} k_2 ε_{eff}}{\\sqrt{\\Delta x \\Delta y}}
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
\\dot{c} = I H_c^{\\alpha}
\\end{aligned}
```

# Fields
 - `I::T`: Calving coefficient
 - `α::T`: Exponent
 - `H_critical::T`: Critical ice thickness above sea level for calving to occur
 - `max_rate::T`: Maximum calving rate

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
\\dot{c} = \\dfrac{f_{ice} max(H_{eff} - H_{max}, 0)}{\\tau_c}
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

Calving law based on the bedrock standard deviation.
```math
\\begin{aligned}
\\dot{c} = \\mathrm{min}\\left(c_{max} f_{scale}, \\dfrac{H_{eff}}{\\tau_c}\\right)
\\end{aligned}
```

# Fields
 - `sd_min::T`: Minimum bedrock standard deviation
 - `sd_max::T`: Maximum bedrock standard deviation
 - `timescale::T`: Calving time scale
 - `max_rate::T`: Maximum calving rate
"""
@kwdef struct BedStddevCalving{T} <: AbstractCalving
    sd_min::T = 50.0            # m
    sd_max::T = 300.0           # m
    timescale::T = 1.0          # yr
    max_rate::T = Inf           # m yr-1
end

############################################################################
# Functions
############################################################################

"""
$(TYPEDSIGNATURES)

Calculate the calving flux based on the ice thickness `H` and the calving law `calving<:AbstractCalving`.
"""
function calving_rate(H, calving::ConstantCalving)
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

function calving_rate(H_eff, tau1, tau2, calving::LipscombCalving)
    
    (; k_τ, w2, max_rate) = calving
    tau_eff = sqrt(max(tau1, 0)^2 + w2 * max(tau2, 0)^2)
    rate = H_eff * k_τ * tau_eff / sqrt(dx * dy)
    return -saturate(rate, 0, max_rate)
end

function calving_rate(H_eff, eps_eff, calving::LevermannCalving)

    (; k2, max_rate) = calving
    calving_ref = max(k2 * eps_eff, 0)
    rate = H_eff * calving_ref / sqrt(dx * dy)
    return -saturate(rate, 0, max_rate)
end

function calving_rate(H, calving::CrawfordCalving)

    (; I, α, H_critical, max_rate) = calving
    H_c = z_srf - z_sl
    if H_c > H_critical && z_bed < z_sl
        c_dt = I * H_c ^ α
    else
        c_dt = 0.0
    end
    return -saturate(c_dt, 0, max_rate)
end

function calving_rate(H, calving::BassisCalving)
    (; C0, α, r, max_rate) = calving
    H_ocn_now = seawater_depth(ρ_seawater__div__ρ_ice, H_ocn)
    timescale = C0 + 0.5 * α * ρ_ice__tim__g * H_eff
    H_max = (1 - r) * timescale / ρ_ice__tim__g +
        sqrt((1 - r)^2 * (timescale / ρ_ice__tim__g)^2 + ρ_seawater__div__ρ_ice * H_ocn ^ 2)

    if H_eff <= H_max
        return 0.0
    else
        return f_ice * max(H_eff - H_max, 0) / timescale
    end
end

function calving_rate(H, calving::BedStddevCalving)
    (; sd_min, sd_max, max_rate, timescale) = calving
    f_scale = (z_bed_stddev - sd_min) / (sd_max - sd_min)
    if f_scale < 0
        f_scale = 0.0
    elseif f_scale > 1
        f_scale = 1.0
    end
    c_dt = min(c_max * f_scale, H_eff / timescale)
    return -saturate(
        c_dt,
        0,
        calving.max_rate,
    )
end

function calving_rate(H::A, c::AbstractCalving) where {A<:AbstractArray}
    c_dt = similar(H)
    calving_rate!(H, c_dt, c)
    return c_dt
end


"""
$(TYPEDSIGNATURES)

Update the calving flux `c_dt` based on the calving law `calving<:AbstractCalving` and the ice thickness `H`.
"""
function calving_rate!(H, c_dt, calving::AbstractCalving)
    map!(h -> calving_rate(h, calving), c_dt, H)
    return nothing
end
###########################################################
# Structs
###########################################################

"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the creep function computation via [`creep`](@ref) and [`creep!`](@ref).
"""
abstract type AbstractCreep end

"""
$(TYPEDSIGNATURES)

Glen-Nye creep function following [glen_creep_1955](@citet) and [nye_distribution_1957](@citet):

```math
\\begin{aligned}
f(\\sigma_e) = E \\, \\sigma_e^{n - 1}
\\end{aligned}
```

The enhancement factor ``E`` accounts for crystal anisotropy, impurity content, or
fabric-induced softening/hardening. ``E > 1`` softens the ice; ``E < 1`` hardens it.

# Fields
 - `n::T=3.0`: Glen-Nye's flow law exponent.
 - `E::T=1.0`: enhancement factor (dimensionless).
"""
@kwdef struct GlenNyeCreep{T} <: AbstractCreep
    n::T = 3.0
    E::T = 1.0
end

"""
$(TYPEDSIGNATURES)

Regularized Glen-Nye creep function following [glen_creep_1955](@citet) and [nye_distribution_1957](@citet):

```math
\\begin{aligned}
f(\\sigma_e) = E \\left( \\sigma_e^{n - 1} + \\sigma_0^{n - 1} \\right)
\\end{aligned}
```

The floor ``\\sigma_0^{n-1}`` prevents the viscosity singularity at zero stress.
In practice ``\\sigma_0 \\ll \\sigma_e`` for realistic stress levels, so the
regularization is only active near ice divides where stresses approach zero.

# Fields
 - `n::T=3.0`: Glen-Nye's flow law exponent.
 - `σ_0::T=1e-6`: regularization stress (``\\mathrm{Pa}``).
 - `E::T=1.0`: enhancement factor (dimensionless).
"""
@kwdef struct RegularizedGlenNyeCreep{T} <: AbstractCreep
    n::T = 3.0
    σ_0::T = 1e-6
    E::T = 1.0
end

"""
$(TYPEDSIGNATURES)

Polynomial creep function following [smith_viscous_1981](@citet):

```math
\\begin{aligned}
f(\\sigma_e) = \\frac{D_0}{\\sigma_0}
\\left( p_0 + p_2 \\left(\\frac{\\sigma_e}{\\sigma_0}\\right)^{\\!2}
           + p_4 \\left(\\frac{\\sigma_e}{\\sigma_0}\\right)^{\\!4} \\right)
\\end{aligned}
```

The polynomial is fit to laboratory data. Unlike the Glen-Nye family, this
parameterization is calibrated to pair with [`SmithMorlandRateFactor`](@ref)
via [`SmithMorlandFlowLaw`](@ref); see that docstring for units and caveats.

!!! warning
    The field names `p1`, `p2`, `p3` in the docstring differ from the actual struct
    fields `p0`, `p2`, `p4` — use the struct fields directly.

# Fields
 - `p0::T=0.3336`: zeroth-order polynomial coefficient.
 - `p2::T=0.3200`: second-order polynomial coefficient.
 - `p4::T=0.02963`: fourth-order polynomial coefficient.
 - `D_0::T=3.169e-8`: strain-rate pre-factor (``\\mathrm{Pa}^{-1}\\,\\mathrm{s}^{-1}``).
 - `σ_0::T=1e5`: reference stress (``\\mathrm{Pa}``).
"""
@kwdef struct SmithMorlandCreep{T} <: AbstractCreep
    p0::T = 0.3336
    p2::T = 0.3200
    p4::T = 0.02963
    D_0::T = 3.169e-8
    σ_0::T = 1e5
end

"""
$(TYPEDSIGNATURES)

Composite creep function following [pettit_ice_2003](@citet), combining a linear
(Newtonian) and a nonlinear (Glen-type) contribution:

```math
\\begin{aligned}
f(\\sigma_e) = A_{\\mathrm{lin}} + A_{\\mathrm{nl}} \\, \\sigma_e^{n - 1}
\\end{aligned}
```

At low stresses the linear term dominates and the viscosity approaches the
Newtonian limit; at high stresses the Glen power-law term takes over.
Reduces to [`GlenNyeCreep`](@ref) when ``A_{\\mathrm{lin}} = 0``.
Suitable for modelling flow near ice divides.

# Fields
 - `n::T=3.0`: Glen-type flow exponent.
 - `A_lin::T=0.0`: linear (Newtonian) coefficient (``\\mathrm{Pa}^{n-1}``).
 - `A_nl::T=1.0`: nonlinear (Glen) coefficient (dimensionless).
"""
@kwdef struct PettitWaddingtonCreep{T} <: AbstractCreep
    n::T = 3.0
    A_lin::T = 0.0
    A_nl::T = 1.0
end

"""
$(TYPEDSIGNATURES)

Composite creep function following [goldsby_superplastic_2001](@citet). Diffusion
creep, grain-boundary sliding (GBS), and basal dislocation slip operate in parallel,
with GBS and basal slip coupled in series:

```math
\\begin{aligned}
\\dot{\\varepsilon}_{\\mathrm{diff}}  &= A_{\\mathrm{diff}} \\, d^{-2} \\, \\sigma_e \\\\
\\dot{\\varepsilon}_{\\mathrm{gbs}}   &= A_{\\mathrm{gbs}}  \\, d^{-p_{\\mathrm{gbs}}} \\, \\sigma_e^{n_{\\mathrm{gbs}}} \\\\
\\dot{\\varepsilon}_{\\mathrm{basal}} &= A_{\\mathrm{basal}} \\, \\sigma_e^{n_{\\mathrm{basal}}} \\\\[4pt]
\\dot{\\varepsilon}_{\\mathrm{tot}}   &= \\dot{\\varepsilon}_{\\mathrm{diff}}
    + \\left(\\dot{\\varepsilon}_{\\mathrm{gbs}}^{-1} + \\dot{\\varepsilon}_{\\mathrm{basal}}^{-1}\\right)^{-1}
\\end{aligned}
```

The creep function returned is ``f(\\sigma_e) = \\dot{\\varepsilon}_{\\mathrm{tot}} / \\sigma_e``.
The pre-factors ``A_{\\mathrm{diff}}``, ``A_{\\mathrm{gbs}}``, ``A_{\\mathrm{basal}}``
carry their own Arrhenius temperature dependence and must be evaluated at the local
temperature before constructing this struct.

# Fields
 - `d::T`: mean grain size (``\\mathrm{m}``).
 - `A_diff::T`: diffusion creep pre-factor (``\\mathrm{Pa}^{-1}\\,\\mathrm{m}^2\\,\\mathrm{s}^{-1}``).
 - `A_gbs::T`: GBS pre-factor (``\\mathrm{Pa}^{-n_{\\mathrm{gbs}}}\\,\\mathrm{m}^{p_{\\mathrm{gbs}}}\\,\\mathrm{s}^{-1}``).
 - `A_basal::T`: basal slip pre-factor (``\\mathrm{Pa}^{-n_{\\mathrm{basal}}}\\,\\mathrm{s}^{-1}``).
 - `p_gbs::T=1.4`: grain-size exponent for GBS.
 - `n_gbs::T=1.8`: stress exponent for GBS.
 - `n_basal::T=2.4`: stress exponent for basal slip.
"""
@kwdef struct GoldsbyKohlstedtCreep{T} <: AbstractCreep
    d::T
    A_diff::T
    A_gbs::T
    A_basal::T
    p_gbs::T = 1.4
    n_gbs::T = 1.8
    n_basal::T = 2.4
end

###########################################################
# Dispatch
###########################################################

"""
$(TYPEDSIGNATURES)

Compute the creep function value based on the effective stress `σ_e` and the creep function parameterization `law<:AbstractCreep`.
"""
function creep(
    σ_e::T,
    law::SmithMorlandCreep,
) where {T<:Real}
    (; p0, p2, p4, D_0, σ_0) = law
    return D_0 / σ_0 * (p0 + p2 * (σ_e / σ_0)^2 + p4 * (σ_e / σ_0)^4)
end

function creep(
    σ_e::T,
    law::GlenNyeCreep,
) where {T<:Real}
    return law.E * σ_e^(law.n - 1)
end

function creep(
    σ_e::T,
    law::RegularizedGlenNyeCreep,
) where {T<:Real}
    (; n, σ_0) = law
    return law.E * (σ_e^(n - 1) + σ_0^(n - 1))
end

function creep(
    σ_e::T,
    law::PettitWaddingtonCreep,
) where {T<:Real}
    (; n, A_lin, A_nl) = law
    return A_lin + A_nl * σ_e^(n - 1)
end

function creep(
    σ_e::T,
    law::GoldsbyKohlstedtCreep,
) where {T<:Real}
    (; d, A_diff, A_gbs, A_basal, p_gbs, n_gbs, n_basal) = law
    ε̇_diff  = A_diff * d^(-2) * σ_e
    ε̇_gbs   = A_gbs * d^(-p_gbs) * σ_e^n_gbs
    ε̇_basal = A_basal * σ_e^n_basal
    ε̇_eff   = ε̇_diff + inv(inv(ε̇_gbs) + inv(ε̇_basal))
    return ε̇_eff / σ_e
end

function creep(
    σ_e::T,
    law::GlenNyeCreep,
) where {T<:Real}
    return law.E * σ_e^(law.n - 1)
end

function creep(
    σ_e::M,
    law,
) where {M<:AbstractArray}
    cf = similar(σ_e)
    creep!(cf, σ_e, law)
    return cf
end

"""
$(TYPEDSIGNATURES)

Same as [`creep`](@ref) but operates in place.
"""
function creep!(cf, σ_e, law)
    map!(s -> creep(s, law), cf, σ_e)
    return nothing
end
###########################################################
# Structs
###########################################################

"""
$(TYPEDSIGNATURES)

An abstract type to multiple dispatch the fabric evolution via [`anisotropy!`](@ref)
and the enhancement factor via [`enhancement_factor`](@ref).
"""
abstract type AbstractAnisotropy end

"""
$(TYPEDSIGNATURES)

Scalar (constant) enhancement factor. The second-order orientation tensor is not
evolved; the enhancement factor is applied uniformly.

# Fields
 - `E::T`: constant enhancement factor (dimensionless).
"""
@kwdef struct EnhancementAnisotropy{T} <: AbstractAnisotropy
    E::T = 1.0
end

"""
$(TYPEDSIGNATURES)

CAFFE (Continuum-mechanical Anisotropic Flow model based on an anisotropic Flow
Enhancement factor) following [placidi_caffe_2010](@citet) and [seddik_caffe_2011](@citet).

The fabric state is carried by the second-order orientation tensor ``\\mathbf{A}``,
defined as the ensemble average of the dyadic products of individual c-axis unit
vectors ``\\mathbf{c}``:

```math
\\mathbf{A} = \\langle \\mathbf{c} \\otimes \\mathbf{c} \\rangle
```

``\\mathbf{A}`` is symmetric, positive semi-definite and satisfies
``\\mathrm{tr}(\\mathbf{A}) = 1``. Its eigenvalues ``(\\lambda_1, \\lambda_2, \\lambda_3)``
characterise the fabric:

| Fabric type | Eigenvalues |
|:---|:---|
| Isotropic | ``\\lambda_i = 1/3`` |
| Single maximum | ``\\lambda_1 \\approx 1, \\; \\lambda_2, \\lambda_3 \\approx 0`` |
| Girdle | ``\\lambda_3 \\approx 0, \\; \\lambda_1 \\approx \\lambda_2 \\approx 1/2`` |

The fabric evolves under the velocity gradient via [`anisotropy!`](@ref).
The flow-law enhancement factor ``E`` for a prescribed loading direction ``\\hat{n}``
is computed by [`enhancement_factor`](@ref).

# Fields
 - `α::T=0.06`: fabric–deformation coupling parameter ``\\alpha \\in [0, 1]``.
   ``\\alpha = 0`` rotates c-axes as rigid bodies (pure spin); ``\\alpha = 1`` gives
   fully affine deformation (isotropic redistribution). Calibrated to ``0.06`` in
   [seddik_caffe_2011](@citet).
 - `E_s::T=10.0`: enhancement factor for a single-maximum fabric perfectly aligned with
   the loading direction (soft limit, ``E_s > 1``).
 - `E_c::T=1.0`: enhancement factor for isotropic fabric (``\\mathrm{tr}(\\mathbf{A}) = 1/3``
   projected onto ``\\hat{n}``). Typically ``E_c = 1``.
"""
@kwdef struct CAFFEAnisotropy{T} <: AbstractAnisotropy
    α::T = 0.06
    E_s::T = 10.0
    E_c::T = 1.0
end

@kwdef struct SimpleCAFFEAnisotropy{T} <: AbstractAnisotropy
    E_min::T = 0.1
    E_max::T = 10.0
    t::T = 8/21 * (E_max - 1) / (1 - E_min)
end

###########################################################
# Dispatch
###########################################################

"""
$(TYPEDSIGNATURES)

Compute the time rate of change of the second-order orientation tensor ``\\mathbf{A}``
and store it in `ȧ_2`. Integrating this rate advances the fabric state.
"""
function anisotropy!(ȧ_2, a_2, D, W, law::CAFFEAnisotropy)
    α = law.α
    trAD = dot(a_2, D)   # tr(A·D) = Frobenius inner product (A symmetric, D symmetric)
    ȧ_2 .= W * a_2 - a_2 * W + (1 - α) .* (D * a_2 + a_2 * D - 2trAD .* a_2)
    return nothing
end

"""
$(TYPEDSIGNATURES)

Compute the CAFFE enhancement factor ``E`` for a given orientation tensor `a_2` and
a unit loading direction `n̂`:

```math
E(\\mathbf{A}, \\hat{n}) = E_c + (E_s - E_c)
\\left(\\frac{3 \\, \\hat{n}^{\\mathrm{T}} \\mathbf{A} \\hat{n} - 1}{2}\\right)^{\\!2}
```

The alignment scalar ``a = \\hat{n}^{\\mathrm{T}} \\mathbf{A} \\hat{n} \\in [0, 1]``
measures the projection of the c-axis cluster onto the loading direction. The factor
``(3a - 1)/2`` is the second Legendre polynomial ``P_2(\\sqrt{a})``, which equals 0 for
isotropic fabric (``a = 1/3``) and 1 for a perfectly aligned single maximum (``a = 1``).

| Fabric state | ``a`` | ``E`` |
|:---|:---|:---|
| Isotropic | ``1/3`` | ``E_c`` |
| Single max, aligned | ``1`` | ``E_s`` |
| Girdle, perpendicular | ``0`` | ``E_c + (E_s - E_c)/4`` |

# Arguments
 - `a_2`: orientation tensor ``\\mathbf{A}`` (3×3 symmetric matrix).
 - `n̂`: unit loading direction vector (length-3).
 - `law`: [`CAFFEAnisotropy`](@ref) parameterisation.
"""
function enhancement_factor(a_2, n̂, law::CAFFEAnisotropy)
    (; E_s, E_c) = law
    a = dot(n̂, a_2 * n̂)
    ψ = (3a - 1) / 2
    return E_c + (E_s - E_c) * ψ^2
end

function enhancement_factor(A, law::SimpleCAFFEAnisotropy)
    (; E_min, E_max, t) = law
    if 0 <= A <= 1
        return E_min + (1 - E_min) * A^t
    elseif 1 < A <= 5/2
        return (4*A^2*(E_max - 1) + 25 - 4 * E_max) / 21
    else
        return E_max
    end
end


"""
$(TYPEDSIGNATURES)

Return the constant enhancement factor; `a_2` and `n̂` are ignored.
"""
function enhancement_factor(_, _, law::EnhancementAnisotropy)
    return law.E
end


"""
$(TYPEDSIGNATURES)

Compute the deformability scalar from the loading direction `n` and the deviatoric stress tensor `t_D`:

```math
\\mathcal{A}(\\hat{n}, \\mathbf{t}_D) = \\frac{5 \\, \\hat{n}^{\\mathrm{T}} \\mathbf{t}_D \\hat{n} - \\mathrm{tr}(\\mathbf{t}_D)}{\\mathrm{tr}(\\mathbf{t}_D)}}
```
"""
function deformability(n, t_D)
    return 5 * square_tangential_invariant(n, t_D) / tr(t_D) ^ 2
end

"""
$(TYPEDSIGNATURES)

Compute the square of the tangential invariant from the loading direction `n` and the deviatoric stress tensor `t_D`:
```math
\\mathcal{A}^2(\\hat{n}, \\mathbf{t}_D) = \\hat{n}^{\\mathrm{T}} \\mathbf{t}_D \\hat{n} - (\\hat{n}^{\\mathrm{T}} \\mathbf{t}_D \\hat{n})^2
```
"""
function square_tangential_invariant(n, t_D)
    return (t_D * n) - (n * t_D * n) ^ 2
end


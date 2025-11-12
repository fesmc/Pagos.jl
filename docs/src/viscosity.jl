#=

# [Flow law](@id flow_law)

The deformation of polycrystalline ice under shear stress ``\tau`` can be described by a flow law that depends on the ice viscosity ``\eta`` and yields the shear (deformation) rate ``\dot{\gamma}``, which can be expressed as:

```math
\begin{aligned}
\dot{\gamma} = \dfrac{1}{\eta(T, p, \tau)} \, \tau
\end{aligned}
```

The viscosity itself depends on temperature (``T``), pressure (``p``), and shear stress (``\tau``). It can be factorized as:

```math
\begin{aligned}
\eta(T, p, \tau) = \frac{1}{2 \, A(T, p) f(|\tau|)}
\end{aligned}
```

where ``A(T, p)`` is the rate factor, and ``f(|\tau|)`` is the creep function. Let's first look at the rate factor.

## [Rate Factor](@id rate_factor)

The rate factor ``A(T, p)`` captures the dependence of ice deformation on temperature and pressure. It typically increases with temperature and decreases with pressure, reflecting the fact that warmer ice deforms more easily, while higher pressure tends to inhibit deformation.

### [Arrhenius Law](@id arrhenius_law)

The rate factor ``A(T, p)`` is typically modeled using an Arrhenius law:

```math
\begin{aligned}
A(T, p) = A_0 \, e^{-\frac{Q}{R T'}}
\end{aligned}
```

where:
- ``A_0`` is the pre-exponential factor,
- ``Q`` is the activation energy,
- ``R`` is the universal gas constant,
- ``T'`` is the temperature relative to the pressure melting point.

For ice, the rate factor is defined in a piecewise manner based on temperature ranges:
```math
\begin{aligned}
A(T, p) =
\begin{cases}
A_1 \, e^{-\frac{Q_1}{R T'}} & \text{for } T' < T^* \\
A_2 \, e^{-\frac{Q_2}{R T'}} & \text{for } T' \geq T^*
\end{cases}
\end{aligned}
```

To define this behaviour in Pagos.jl, use [`ArrheniusRateFactor`](@ref):
=#

using Pagos, CairoMakie
arrhenius_rate_factor = ArrheniusRateFactor()
Tprime = range(-50, stop = 10, length = 121) .+ 273.15
A = map(x -> get_rate_factor(x, arrhenius_rate_factor), Tprime)
fig = plot_rate_factor(Tprime .- 273.15, A)

#=

This matches Fig. 4.5 of [greve_dynamics_2009](@citet)! In [`AbstractRateFactor`](@ref), we show other options, as well as how to implement your own rate factor models.

## [Creep function](@id creep_function)

The creep function ``f(|\tau|)`` describes how the viscosity changes with applied shear stress. A commonly used creep function is the Glen-Nye flow law, which is defined as:

```math
\begin{aligned}
f(|\tau|) = |\tau|^{n - 1}
\end{aligned}
```

where ``n`` is the creep exponent, typically around 3 for ice. This means that the viscosity decreases with increasing shear stress, leading to non-linear deformation behavior.

To ease the implementation regardless of the coordinate system, the effective stress ``\sigma_e`` is often used. It is defined as the square root of the second invariant of the stress tensor and simplifies the creep function to:

```math
\begin{aligned}
f(\sigma_e) = \sigma_e^{n - 1}
\end{aligned}
```

Thus, the viscosity equation can be rewritten as:
```math
\begin{aligned}
\eta(T', \sigma_e) = \frac{1}{2 \, A(T') \, \sigma_e^{n - 1}}
\end{aligned}
```

However, this is diverges for ``\sigma_e \to 0``. To avoid this, a regularization parameter ``\sigma_0`` is often introduced:

```math
\begin{aligned}
\eta(T', \sigma_e) = \frac{1}{2 \, A(T') \, (\sigma_e + \sigma_0)^{n - 1}}
\end{aligned}
```

This is known as the Regularized Glen-Nye flow law, and is implemented in Pagos.jl as [`RegularizedGlenNyeViscosity`](@ref):
=#

glennye_flowlaw = RegularizedGlenNyeViscosity()
A = 1e-24
σ_e = 100f3
η_ice = get_viscosity(A, σ_e, glennye_flowlaw)

#=
This exemplifies the computation for a single value of effective stress and rate factor. Let's now look at how the viscosity varies with effective stress for different temperatures and reproduce Fig. 4.6 of [greve_dynamics_2009](@citet).
=#

T = (0, -10, -20)
A = map(x -> get_rate_factor(x + 273.15, arrhenius_rate_factor), T)
σ_e = range(0, stop = 100, step = 0.1) .* 1f3
η_ice = [map(x -> get_viscosity(AA, x, glennye_flowlaw), σ_e) for AA in A]
fig = plot_ice_viscosity(η_ice, σ_e, T)

#=
In [`AbstractViscosity`](@ref), we show other options, as well as how to implement your own viscosity models.

## [Temperature relative to pressure melting point](@id melting_point)

For some computations, it is necessary to determine the temperature relative to the pressure melting point, ``T'``. The pressure melting point decreases with increasing pressure, and can be approximated using a linear relation:
=#

lpmp = LinearPressureMeltingPoint()
p = range(0, stop = 50f6, length = 1000)
Tm = map(x -> get_melting_point(x, lpmp), p)
fig_tm = plot_melting_point(p, Tm)

#=
In [`AbstractPressureMeltingPoint`](@ref), we show other options, as well as how to implement your own pressure melting point models.
=#
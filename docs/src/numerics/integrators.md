```@meta
EditURL = "integrators.jl"
```

# [Time integrators](@id integrators)

Pagos ships a small, self-contained set of **explicit** time integrators. They are
deliberately minimal — just enough to step a component forward in time — with no
dependency on DifferentialEquations.jl. All of them:

- are **matrix-free** and run unchanged on CPU and GPU (the state arrays decide the backend),
- are **allocation-free** in the hot loop (all scratch is pre-allocated at construction),
- are **type-stable** (verified with `@inferred`),
- share one interface, so you can swap methods without touching the rest of your code.

This page is a decision guide: *"I have this kind of problem → which method should I pick?"*
For the full list of types and fields, see the [`AbstractIntegrationMethod`](@ref) API.

## The common interface

You provide an in-place right-hand side `f(du, x, p, t)` (a top-level `function` or a
callable struct — *not* a capturing closure), wrap it together with the state `x`,
parameters `p`, and a method in an [`Integrator`](@ref), and advance it with
[`step!`](@ref):

````@example integrators
using Pagos

# RHS of the form f(du, x, p, t). Example: exponential decay dx/dt = -x.
decay!(du, x, p, t) = (@. du = -x)

x     = [1.0, 2.0, 3.0]
integ = Integrator(decay!, copy(x), nothing, Tsitouras54(x; dt = 0.1))
step!(integ, 1.0)          # advance by a macro-step Δt_sync = 1.0
integ.x                    # ≈ x .* exp(-1)
````

`step!(integ, Δt_sync)` advances the state by a **macro-step** `Δt_sync`, taking as many
internal sub-steps as the method needs. The right-hand side is treated as *autonomous*
over the window (boundary conditions frozen for the duration of `Δt_sync`), which matches
the operator-splitting coupling used in `step!(::Simulation, …)`. You may freely mutate
`x` and `p` between macro-steps.

Because each physics component owns its own `time_stepper`, you can pick a **different
integrator per component** — e.g. an SSP method for thickness transport and a stabilized
method for the thermal step.

## Which method should I choose?

| You have… | Use | What controls the step |
|---|---|---|
| A smooth, non-stiff problem (the default) | [`Tsitouras54`](@ref) (Tsit5) | accuracy: `atol`, `rtol` |
| A cheap RHS, or you only need crude accuracy | [`BogackiShampine32`](@ref) (BS32) | accuracy: `atol`, `rtol` |
| A need for a fixed, reproducible step | [`RungeKutta4`](@ref) / [`Euler`](@ref) | fixed `dt` |
| A diffusion-dominated / stiff problem (SIA thickness, vertical heat) | [`RKL2`](@ref) | `dt` + `spectral_radius` |
| Advection of a non-negative field with sharp fronts (thickness) | [`SSPRK43`](@ref) / [`SSPRK33`](@ref) | CFL: `dt_fe` |
| A genuinely stiff problem needing unconditional stability | *not in this toolbox* (use an implicit / IMEX solver) | — |

The two questions that decide almost everything are: **is the problem stiff?** and
**what must the step size respect — accuracy, stability, or positivity?**

### "I just want something that works"

Use [`Tsitouras54`](@ref). It is an efficient 5th-order embedded pair with automatic,
error-controlled step sizing — a strong default for smooth, non-stiff right-hand sides.
Tighten or loosen `atol`/`rtol` to trade accuracy for speed:

````@example integrators
tsit = Tsitouras54(x; dt = 0.1, atol = 1e-8, rtol = 1e-8)   # initial dt; tolerances drive adaptivity
````

### "My RHS is cheap, or I only need low accuracy"

[`BogackiShampine32`](@ref) is a 3rd-order embedded pair — fewer stages per step than
Tsit5, so it wins when evaluations are cheap or the tolerance is loose.

### "I want a fixed step (debugging, reproducibility, teaching)"

[`RungeKutta4`](@ref) (4th order) or [`Euler`](@ref) (1st order) take a fixed `dt` with no
error control. Predictable and simple, but you are responsible for choosing a stable `dt`.

````@example integrators
rk4 = RungeKutta4(x; dt = 1e-2)
````

### "My problem is diffusion-dominated / stiff"

The shallow-ice (SIA) thickness equation is a *nonlinear diffusion*, and vertical heat
conduction is linear diffusion — both stiff. An accuracy-controlled explicit pair
(Tsit5/BS32) will not blow up, but it gets throttled to the tiny *stability* step
``\Delta t \sim \Delta x^2 / D`` and crawls.

[`RKL2`](@ref) (Runge–Kutta–Legendre super-time-stepping) is built for exactly this. It
stays **explicit and matrix-free** but covers a superstep with ``s`` cheap stages whose
stability region stretches along the negative real axis
(``\Delta t_\mathrm{stable} \approx \Delta t_\mathrm{expl}\,(s^2+s-2)/4``), so it takes far
larger steps than Tsit5 on diffusive problems. Its adaptivity is in the *stage count*
``s``, not the step size — `dt` is a fixed target superstep.

RKL2 needs the spectral radius ``\rho = \max|\lambda(\partial f/\partial x)|``. You can let
it estimate ``\rho`` automatically each superstep (matrix-free power iteration), or pass an
analytic bound if you have one (cheapest and safest):

````@example integrators
rkl2_auto  = RKL2(x; dt = 0.1)                       # ρ estimated automatically
rkl2_fixed = RKL2(x; dt = 0.1, spectral_radius = 4.0) # ρ supplied (e.g. 2D(1/Δx²+1/Δy²) for diffusivity D)
````

!!! tip "Updating ρ"
    The spectral radius depends on the state. In auto mode it is refreshed every
    superstep; if you supply it, update `method.spectral_radius` whenever the state or
    parameters change. Underestimating ``\rho`` is the *unsafe* direction — when in doubt,
    bound it generously (the auto estimate applies a safety factor).

!!! warning "Super-time-stepping buys stability, not accuracy"
    RKL2 lets you take steps far larger than the explicit stability limit *without blowing
    up*, but it is only 2nd order. Size the superstep to resolve the **dynamical timescale**
    you care about; do not expect accuracy from a single superstep that spans many
    e-folding times.

### "I'm advecting a non-negative field with sharp fronts"

Ice thickness must stay ``\ge 0`` and develops sharp margins. Here the priority is **positivity / monotonicity**, not high order. Strong-stability-
preserving methods preserve any bound that forward Euler satisfies, as long as the step
obeys a CFL limit.

[`SSPRK43`](@ref) (4-stage, 3rd order, SSP coefficient ``C = 2``) is the recommended
choice; [`SSPRK33`](@ref) (3-stage, ``C = 1``) is the classic alternative. You supply the
**forward-Euler-stable (CFL-limited) step** `dt_fe` from the current wave speeds, and the
method safely takes steps of size ``C \cdot dt_\mathrm{fe}``:

````@example integrators
cfl, dx, max_speed = 0.7, 16e3, 1e3      # example: ν·Δx / max|u|
dt_fe = cfl * dx / max_speed
sspp  = SSPRK43(x; dt_fe = dt_fe)        # steps at 2·dt_fe, preserving positivity
ssp_coefficient(sspp)                    # C = 2
````

!!! tip "Updating dt_fe"
    Like RKL2's `spectral_radius`, `dt_fe` is state-dependent — update `method.dt_fe` each
    macro-step from the current velocities, e.g. `dt_fe = cfl·Δx / max|u|`.

### "It's genuinely stiff and I need unconditional stability"

That is outside this explicit toolbox. Unconditionally stable integration needs an
**implicit** or **IMEX** scheme (a linear/nonlinear solve per step). For the stiff thermal
step, the standard remedy is a **vertical-implicit** discretization (a tridiagonal solve
per column), handled by the thermodynamics solver rather than a generic ODE integrator.

## Can I control accuracy in RKL2 and the SSP methods?

[`Tsitouras54`](@ref) and [`BogackiShampine32`](@ref) adapt their step from an embedded
*error* estimate. [`RKL2`](@ref), [`SSPRK33`](@ref) and [`SSPRK43`](@ref) do **not** — and
that is deliberate: for these families the step is normally limited by *stability/CFL*, not
accuracy, so a temporal error controller would usually just sit at the stability cap.

**The knob you already have is the step size**, since all three are fixed order:

- [`RKL2`](@ref) is 2nd order → halving `dt` (the superstep) cuts the error ``\sim 4\times``.
- [`SSPRK33`](@ref)/[`SSPRK43`](@ref) are 3rd order → halving `dt_fe` (or the CFL number)
  cuts the error ``\sim 8\times``.

This is usually the *right* control here: choose the step to resolve the dynamical timescale
(RKL2) or to sit at the CFL limit (SSP), and that choice sets the accuracy. Note two
constraints before reaching for anything fancier:

- For the SSP methods the step is capped at ``C\cdot dt_\mathrm{fe}`` by **positivity**, so
  an accuracy controller could only ever *shrink* below that cap — growing past it would
  forfeit the SSP guarantee. Near sharp fronts the solution is non-smooth anyway, so
  temporal accuracy is limited by the *spatial* scheme.
- [`RKL2`](@ref) handles stability through its *stage count*, not the step, but it carries
  **no embedded estimate** to drive accuracy-based adaptivity.

If you do need *automatic* error control, two routes (neither built in yet):

| Route | Applies to | Cost / notes |
|---|---|---|
| **Step doubling** (one step `h` vs two of `h/2`; difference ``\propto h^{p+1}``) | all three | universal, no special coefficients, ``\sim 3\times`` cost; keep the stability/CFL cap as a hard upper bound, `dt_new = min(error_based, C·dt_fe)` |
| **Embedded pair** (lower-order companion reusing the stages) | [`SSPRK33`](@ref)/[`SSPRK43`](@ref) | efficient, BS32/Tsit5-style control capped at ``C\cdot dt_\mathrm{fe}``; published embedded weights exist (Conde, Fekete & Shadid, 2018) |
| **Switch to ROCK2/ROCK4** | replaces [`RKL2`](@ref) | the stabilized-explicit family that ships with a built-in error estimator |

In short: for most ice-sheet use, control accuracy with the **step size** (and spatial
resolution, plus the CFL number for SSP); add an embedded SSPRK43 controller or ROCK2 only
in a regime where accuracy genuinely binds before stability.

## Cheat sheet

- **Smooth & non-stiff** → [`Tsitouras54`](@ref) (default), [`BogackiShampine32`](@ref) if cheap/crude.
- **Fixed step** → [`RungeKutta4`](@ref) / [`Euler`](@ref).
- **Stiff diffusion (explicit)** → [`RKL2`](@ref), with `spectral_radius` auto or analytic.
- **Positive advection with fronts** → [`SSPRK43`](@ref) / [`SSPRK33`](@ref), with CFL `dt_fe`.
- **Stiff, needs implicit** → not here; use a vertical-implicit / IMEX solver.

A closing reminder on performance: these integrators are bandwidth-bound and their cost is
dominated by your right-hand side `f`. The biggest performance lever is writing `f` well
(e.g. as KernelAbstractions kernels for CPU-threading + GPU + SIMD), not tuning the
integrator itself.


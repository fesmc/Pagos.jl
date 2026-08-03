# ---------------------------------------------------------------------------
# Time integrators (`src/utils/integrators.jl`).
#
# `step!` advances the state over a macro-step by taking as many internal sub-steps as the
# method needs, so what is timed is *one macro-step of fixed size* — same `Δt_sync`, same
# `dt`, hence the same number of stage evaluations for every method and every sample. The
# right-hand side is deliberately cheap and shared, so the entries compare scheme overhead
# (stage count, buffer traffic) rather than RHS cost.
#
# Only fixed-step methods are tracked. The adaptive ones (`BogackiShampine32`,
# `Tsitouras54`) write the accepted step size back into `m.dt` between calls, so sample n
# would start from a different `dt` than sample 1 and the entry would measure a moving
# target. `RKL2` is included with an explicit `spectral_radius`, which pins its stage count
# and skips the power iteration for the same reason.
# ---------------------------------------------------------------------------

# Non-capturing, top-level, `f(du, x, p, t)` — the form `Integrator` documents. A closure
# over the work array would be captured by value into every stage call.
function _diffusion_rhs!(du, u, p, t)
    ∂x!(p.work, u, p.dx)
    ∂x!(du, p.work, p.dx)
    return nothing
end

"""
    integrators_suite(fx)

One macro-step of `step!` per fixed-step scheme, over a 2D field with a diffusion-shaped
right-hand side.
"""
function integrators_suite(fx)
    (; backend, T, nx, ny, dx) = fx
    g = BenchmarkGroup()

    host = T[sinpi(2i / nx) * cospi(2j / ny) for i in 1:nx, j in 1:ny]
    x0   = to_backend(backend, host)
    p    = (; work = to_backend(backend, zero(host)), dx = T(dx))

    dt = T(1.0e-3)
    # One macro-step = four sub-steps for every method, so the entries stay comparable.
    dt_sync = 4dt

    # `setup` restores the state before each sample: `step!` mutates `x` in place, and an
    # unbounded diffusion march would drift the field (and eventually its magnitude) away
    # from what sample 1 saw.
    for (name, make) in (
            ("Euler",    x -> Euler(x; dt = dt)),
            ("RK4",      x -> RungeKutta4(x; dt = dt)),
            ("SSPRK33",  x -> SSPRK33(x; dt_fe = dt)),
            ("SSPRK43",  x -> SSPRK43(x; dt_fe = dt)),
            ("RKL2",     x -> RKL2(x; dt = dt, spectral_radius = 1 / dt)),
        )
        g[name] = @benchmarkable begin
            step!(integ, $dt_sync)
            sync!($backend)
        end setup = (
            integ = Integrator(_diffusion_rhs!, copy($x0), $p, $make(copy($x0)))
        ) evals = 1
    end

    return g
end

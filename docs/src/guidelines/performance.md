# Performance

## Asynchronous kernel launches

Pagos operators follow a three-layer pattern:

1. **struct-level wrapper** (e.g. `deviatoric_stress!(mech, mat)`) — unpacks state fields;
2. **array-level launcher** (e.g. `deviatoric_stress!(sxx, ..., eyz)`) — picks the
   backend and launches a KernelAbstractions kernel;
3. **device-side code** — the `@kernel` body and, for extensible operators like
   `strainrate!`, per-element methods dispatching on the momentum balance.

Launchers return after *launching* their kernel, without waiting for it to complete.
This is safe because all kernels, broadcasts and `copyto!`s issued from the same Julia
task are **stream-ordered**: each operation sees the writes of every operation launched
before it. It is also what makes tight solver loops (e.g. the pseudo-transient solver)
fast on GPU — the host enqueues work and runs ahead instead of paying a host–device
round trip per kernel. On the CPU backend, kernel calls complete before returning, so
nothing changes there.

The host must synchronize only when it *reads* device data:

- **Scalar reductions and copies to host** (`maximum(abs, u)`, `Array(u)`, `u[i, j]`)
  synchronize implicitly — no action needed.
- **Timing/benchmarking** must call `KernelAbstractions.synchronize(backend)` inside
  the timed region, otherwise only the launch overhead is measured (see
  `benchmark/stress.jl` for the pattern).
- **Reading from a different Julia task** (a different stream) requires an explicit
  `KernelAbstractions.synchronize(backend)` before handing over the data.

For **external users extending Pagos**: extension points are per-element methods
(layer 3), which contain pure stencil math and no synchronization logic — new physics
inherits CPU/GPU support and the async behaviour automatically. Only when authoring an
entirely *new* operator (a new launcher) does the convention apply: do not synchronize
inside the launcher; leave it to the caller. Adding a `synchronize` anyway is never
incorrect, it merely reintroduces a stall.

## Systematical benchmark

In future, the documentation will include some benchmarks. Every PR that modifies existing routines needs to prove that the code is not slowed down or, if so, for a good reason.

## Using `map!`

The in-place `map!` applies functions to arrays without allocation. It combines well with different array types, e.g., `SubArray`, `CuArray`, etc. It might be less performant for specific cases but is much more general and should therefore be preferred. Simple benchmarks showed that applying `map!` over the full array is faster than applying it only to a view of the array. Therefore, if some operation is not defined outside of a mask, it is better to apply `map!` to the full array and use a conditional inside the function.
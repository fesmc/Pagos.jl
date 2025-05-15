
# Pagos.jl - Πάγος

<!-- [![Build Status](https://github.com/fesmc/Pagos.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fesmc/Pagos.jl/actions/workflows/CI.yml?query=branch%3Amain) -->

## Aim

Pagos.jl is an ice-sheet model written in pure Julia. It is designed to be accessible,
modular and performant. We believe that these are the key aspects of a research model
since it allows a widespread use, participative development, a simplified permutation
of parameters and laws and a reduced time-to-first-paper. To this end, we expect that as much care will be put into designing the numerics and internal procedures used in the model as into its API and usability.

If you want to go to multi-GPU computation and unlock very high resolutions, we recommend
to look at [FastIce.jl](https://github.com/PTsolvers/FastIce.jl). If you want to apply
autodifferentiation (AD), we recommend to look at [dJuice](https://github.com/DJ4Earth/dJUICE.jl).
Although these two projects comply with specific requirements, they might lack
the representation of important processes that are included in Pagos.jl. In the future,
Pagos.jl will be compatible with AD (soon after Enzyme.jl and DifferentiationInterface.jl
converge to stable releases).

### Accesibility

Accessibility is the absolute priority of Pagos. This means that the code should be:

- **Legible**
  - Use of long and unambiguous names for functions, variables and structs.
  - Commands should be broken down into steps if it does not imply a performance loss.

- **Documented**
  - Each function or struct comes with extensive docstrings.
  - Reference to papers should be clearly provided using DocumenterCitations.jl.
  - Formulas should be formatted in Latex.
  - Each use case should be illustrated through at least one example.
  - Anything that is not easily understandable or documented should be
explained in greater detail.

- **Open source and open development**
  - The code is open source (distributed under X).
  - The code is easy to install by running `]add Pagos` in the Julia REPL.
  - The development is structured around GitHub issues, which allows people to easily
  request features, join the development process and track past changes.

We emphasise that while modularity and performance can sometimes lead to trade-off
decisions, accessibility can and should always be complied with. In particular, we
provide a guideline for variable names that are recurrent throughout the code. We
invite any developer to thoroughly read and comply with this list, which can be found
in `docs/src/guidelines/naming.md`, as well as with the [julia style guidelines](https://docs.julialang.org/en/v1/manual/style-guide/).

### Modularity

Ice-sheet modelling is subject to a lot of design choices. For instances, various basal
friction laws exist and should be interchangeably usable. In essence, this means that
the user should be able to do:

```julia
using Pagos
rcf = RegularizedCoulombFriction(params...)
ais = IceSheet(args..., friction_law = rcf)
```

and alternatively:

```julia
plf = PowerLawFriction(params...)
ais = IceSheet(args..., friction_law = plf)
```

Modularity also means that the code is leightweight: any module that is not necessary
can be omitted, which leads to faster download and precompilation. The modular view of the code should be held up to date in `docs/src/structure.drawio.png`,
which can be edited via draw.io (including its VS code extension).

### Performance

We distinguish various levels of performance optimization:

- serial performance
- parallel performance
- special architectures

We here focus on serial performance. The reason for this is that 0 memory allocation
and an efficient code on a single CPU should be the baseline for everything else.
It should already be possible to achieve great performance with this, which can
subsequently be easily parallelised through the use of `@threads`. Performance on special
architectures such as GPUs and TPUs are envisioned for the future but not targeted for
the first major release (v1.0).

Performance is achieved through some important techniques that constitute the bread
and butter of the code since they are not standard in ice-sheet modelling yet:

- Pseudo-transient method
- High-order adaptive time stepping (Tsit5)
- FFT-based isostasy
- Wavelet transform?
- FFT-based Picard iteration?

We also target balanced choices concerning the physics to achieve high performance. This means:

- DIVA instead of SIA/SSA vs. Stokes.
- ML techniques can be plugged in easily due to the high modularity of the model.

### Flexibility

- Variables are defined in a general way based on the actual physics of the system, even when it is not necessary, e.g., due to the use of a simplifying approximation. This allows the model to be agnostic to the method producing the variable and facilitates future development.
- Wide range of simple subglacial hydrology models.
- Wide range of fracture models.

### Note for devs

Any new functionality should be described, tested and documented. If not, the pull request
(PR) will not be prioritized.

# Structure

![](/docs/src/structure.drawio.png)

# Benchmark tests

## 1. Slab test

As shown by Robinson et al. (2022).

Test 1 should give a uniform x-velocity of `ux = 8.93` m/yr.
Test 2 should give a uniform x-velocity of `ux = 149` m/yr.

```julia
using Pkg
Pkg.activate(".")
include("test/numerics/slab.jl")
```

## 2. Stream test

As defined by Schoof (2006), and reproduced with different parameter values
by Lipscomb et al. (2019) and Berends et al. (2022).

```julia
using Pkg
Pkg.activate(".")
include("test/numerics/stream.jl")
```

This test will produce a figure showing the analytical stream solution
for the x-velocity and the numerical result for various resolutions.

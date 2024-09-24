
# Pagos.jl - Πάγος

<!-- [![Build Status](https://github.com/JanJereczek/Pagos.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JanJereczek/Pagos.jl/actions/workflows/CI.yml?query=branch%3Amain) -->

# Aim

Pagos.jl is an ice-sheet model written in pure Julia. It is designed to be accessible,
modular and performant. We believe that these are the key aspects of a research model
since it allows a widespread use, a participative development, a simplified permutation
of parameters and laws and a reduced time-to-first-paper.

If you want to go to multi-GPU computation and unlock very high resolutions, we recommend
to look at [FastIce.jl](https://github.com/PTsolvers/FastIce.jl). If you want to apply
autodifferentiation (AD), we recommend to look at [dJuice](https://github.com/DJ4Earth/dJUICE.jl).
Although these two projects comply with specific requirements, they might lack
the representation of important processes that are included in Pagos.jl. In the future,
Pagos.jl will be compatible with AD (soon after Enzyme.jl and DifferentiationInterface.jl
converge to stable releases).

We coin the word interplex, which is a contraction of "intermediate complexity".
An interplex model emphatically goes as high as possible on the scale of complexity
until a minor improvement in accuracy comes at the expense of a significant jump
in complexity. Pagos is designed along this central line to avoid the misrepresentation
of any process while still being able to run relatively high resolutions on paleo time
scales.

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
invite any developper to thoroughly read and comply with this list, which can be found
in `docs/src/naming_guidelines.md`, as well as with the [julia style guidelines](https://docs.julialang.org/en/v1/manual/style-guide/).

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
can be omitted, which leads to faster download and precompilation. For instance, to use
Pagos with a simple quadratic nonlocal melting parametrization for ice shelves:

```julia
qnl = QuadraticNonlocalMelt(params...)
ais = IceSheet(args..., melt = qnl)
```

But more elaborate approaches should be usable just as simply:

```julia
using Laddie
cavity2D = Cavity2D(params...)
ais = IceSheet(args..., melt = cavity2D)
```

The modular view of the code should be held up to date in `docs/src/structure.drawio.png`,
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
architectures such as GPUs and TPUs are wishful for the future but not targetted for
the first major release (v1.0).

### Interplexity

Interplexity is achieved by using:
- DIVA instead of SIA/SSA vs. Stokes,
- FastIsostasy instead of ELRA/ELVA vs. 3D GIA model,
- Laddie instead of QNL/PICO vs. cavity resolving ocean model,
- ML techniques can be plugged in easily due to the high modularity of the model,
- Wide range of simple subglacial hydrology models,
- Wide range of fracture models.

### Numerical bread and butter

Performance is achieved through some important techniques that constitute the bread
and butter of the code since they are not standard in ice-sheet modelling yet:
- Pseudo-transient method
- High-order adaptive time stepping (Tsit5)
- FFT-based isostasy
- Wavelet transform?
- FFT-based Picard iteration?

### Note for devs

Any new functionality should be described, tested and documented. If not, the pull request
(PR) will not be prioritized.

# Structure

![](/docs/src/structure.drawio.png)

# Benchmark tests

## 1. Slab test

As seen in Robinson et al. (2022). 

Test 1 should give a uniform x-velocity of `ux = 8.93` m/yr. 
Test 2 should give a uniform x-velocity of `ux = 149` m/yr. 

```julia
using Pkg
Pkg.activate(".")
include("test/slab.jl")
```

## 2. Stream test

As seen in Schoof (2006), and reproduced with different parameter values
in Lipscomb et al. (2019) and Berends et al. (2022). 

```julia
using Pkg
Pkg.activate(".")
include("test/stream.jl")
```

This test will produce a figure showing the analytical stream solution 
for the x-velocity and the numerical result for various resolutions. 
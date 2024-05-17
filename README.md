
# Pagos.jl - Πάγος

[![Build Status](https://github.com/JanJereczek/Pagos.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JanJereczek/Pagos.jl/actions/workflows/CI.yml?query=branch%3Amain)

A pure Julia ice-sheet model.

# Benchmark tests

## 1. Slab test

As seen in Robinson et al. (2022). 

Test 1 should give a uniform x-velocity of `ux = 8.93` m/yr. 
Test 2 should give a uniform x-velocity of `ux = 149` m/yr. 

```
using Pkg; Pkg.activate(".")

include("test/slab.jl")
```

## 2. Stream test

As seen in Schoof (2006), and reproduced with different parameter values
in Lipscomb et al. (2019) and Berends et al. (2022). 

```
using Pkg; Pkg.activate(".")

include("test/stream.jl")
```

This test will produce a figure showing the analytical stream solution 
for the x-velocity and the numerical result for various resolutions. 
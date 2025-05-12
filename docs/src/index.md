# Getting started

Pagos.jl is an ice-sheet model written in pure Julia. It is designed to be accessible,
modular and performant. To install it, please run:

```
]add https://github.com/JanJereczek/Pagos.jl
```

!!! warning 
Pagos is currently work in progress and is therefore not fully functional and subject to
major changes in near future.

## Modular architecture

The modular architecture of Pagos.jl allows an easy coupling to other models, as well as
lightweight submodules that allow users to merely download the part of the code required
for their application.
"""
$(TYPEDSIGNATURES)

Top-level run configuration wrapping an `IceSheet` with IO and timing.
Separating `Simulation` from `IceSheet` lets restarts swap the output target
without reconstructing the physical state:

```julia
sim  = Simulation(ais)
sim2 = Simulation(ais; io = OutputWriter("restart.jld2"))
```
"""
struct Simulation{IS,IO,TM}
    ice_sheet::IS   # <: IceSheet
    io::IO          # <: AbstractOutputWriter
    timer::TM       # <: AbstractTimer
end

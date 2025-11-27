#=

# [Calving](@id calving)

Calving laws describe the process by which icebergs break off from an ice-sheet margin. This process is influenced by various factors, including ice thickness, temperature, and stress conditions at the ice front. In Pagos.jl, calving laws can be implemented by defining subtypes of the abstract type `AbstractCalving` and providing a method to calculate the calving flux based on ice thickness.

## [No Calving](@id no_calving)

The simplest calving law is the `NoCalving` law, which assumes that no calving occurs. This is represented by the `NoCalving` struct, which is a subtype of `AbstractCalving`. The calving flux in this case is always zero.
=#

no_calving = NoCalving()

#=

## [Threshold Calving](@id threshold_calving)

=#

threshold_calving = ThresholdCalving()
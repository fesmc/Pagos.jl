"""
$(TYPEDSIGNATURES)

Top-level struct representing an ice sheet.

Each field is a component struct that bundles its own state, grid, time stepper,
and physical parameterizations. Components that share the ice sheet's spatial
discretization hold a `CommonGrid` sentinel instead of a separate grid allocation.

Typical usage via a convenience constructor:

```julia
ais = antarctica(;)         # returns an IceSheet
sim = Simulation(ais)
step!(sim, Δt_sync)
```
"""
struct IceSheet{TP, TR, DN, TD, MT, C, CG, BE, PR}
    topography::TP      # <: Topography
    tracers::TR         # <: Tracers
    dynamics::DN        # <: Dynamics
    thermodynamics::TD  # <: Thermodynamics
    material::MT        # <: Material
    constants::C        # <: Constants
    common_grid::CG     # <: AbstractGrid
    backend::BE
    projection::PR
end

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
struct IceSheet{TP,DN,TD,MT,C,CG,BE,PR}
    topography::TP      # <: Topography
    dynamics::DN        # <: Dynamics
    thermodynamics::TD  # <: Thermodynamics
    material::MT        # <: Material
    constants::C        # <: Constants
    common_grid::CG     # <: AbstractGrid
    backend::BE
    projection::PR
end

"""
$(TYPEDSIGNATURES)

# Fields
 - `state`: the topographic state variables
 - `grid`: the grid on which the topography is defined — `Grid`, or `CommonGrid` to use
   the `IceSheet`'s own grid instead of a separate allocation
 - `time_stepper`: `Tsit54()` (or similar), or `CommonTimeStepper()` to reuse the
   `IceSheet`'s stepper and avoid a redundant pre-allocation
 - `mass_balance`: the mass balance model for the topography component
 - `solver`: the mass balance solver for the topography component
 - `calving`: the calving model for the topography component
 - `sigma_transform`: the sigma transform for the topography component
"""
struct Topography{S,G,TS,MB,SL,C,ST}
    state::S            # <: TopographicState
    grid::G             # <: Grid or CommonGrid
    time_stepper::TS    # <: Tsit54() or CommonTimeStepper()
    mass_balance::MB    # <: AbstractMassBalance
    solver::SL          # <: AbstractMassBalanceSolver
    calving::C          # <: AbstractCalving
    sigma_transform::ST # <: AbstractSigmaTransform
end

# @dev TODO: this should move into Ichnos.
"""
$(TYPEDSIGNATURES)

# Fields
 - `state`: the tracer state variables
 - `grid`: the grid on which the tracers are defined
 - `time_stepper`: the time stepper for the tracers component
 - `age`: the age model for the tracers component
 - `isotopes`: the isotopic composition, a `NamedTuple` (`:oxygen_18`, `:deuterium`, ...)
 - `frame`: the reference frame; only `Eulerian` is supported so far
"""
struct Tracers{S,G,TS,A,I,FR}
    state::S            # <: TracersState
    grid::G             # <: Grid or CommonGrid
    time_stepper::TS    # <: Tsit54() or CommonTimeStepper()
    age::A              # <: AbstractAgeModel
    isotopes::I
    frame::FR           # <: AbstractFrame
end

"""
$(TYPEDSIGNATURES)

# Fields
 - `state`: the dynamics state variables
 - `grid`: the grid on which the dynamics are defined
 - `time_stepper`: the time stepper for the dynamics component
 - `momentum`: the momentum balance model for the dynamics component
 - `solver`: the momentum solver for the dynamics component
 - `effective_pressure`: the effective pressure model for the dynamics component
 - `friction`: the basal friction model for the dynamics component
"""
struct Mechanics{S,G,TS,M,SL,EP,F}
    state::S            # <: MechanicState
    grid::G             # <: Grid or CommonGrid
    time_stepper::TS    # <: Tsit54() or CommonTimeStepper()
    momentum::M         # <: AbstractMomentumBalance
    solver::SL          # <: AbstractMomentumSolver
    effective_pressure::EP  # <: AbstractEffectivePressure
    friction::F         # <: BasalFriction
end

"""
$(TYPEDSIGNATURES)

"""
struct Thermodynamics{S,G,TS,E,SL}
    state::S            # <: ThermodynamicsState
    grid::G             # <: Grid or CommonGrid
    time_stepper::TS    # <: Tsit54() or CommonTimeStepper()
    energy::E           # <: AbstractEnergyBalance
    solver::SL          # <: AbstractEnergySolver
end

"""
$(TYPEDSIGNATURES)

"""
struct Material{S,G,TS,A,FL,PMP}
    state::S            # <: MaterialState
    grid::G             # <: Grid or CommonGrid
    time_stepper::TS    # <: Tsit54() or CommonTimeStepper()
    anistropy::A        # <: AbstractAnisotropy
    flow_law::FL        # <: AbstractFlowLaw
    pressure_melting_point::PMP # <: AbstractPressureMeltingPoint
end
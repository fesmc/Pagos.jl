#=
# Kryonomos

Includes:
- Ice sheet
- Surface mass balance
- Basal mass balance from ocean
- Geothermal heat flux
- GIA
- Hydrology
- Barystatic sea level
- Landscape evolution
=#

tspan = (0f0, 1f4)
dt_coupling = 10f0
ice = Pagos.Simulation(..., model)
smb = Chios.SurfaceMassBalance(...)
ioi = Okeanos.IceOceanInterface(...)        # or Pontos
ghf = Tartaros.GeothermalHeatFlux(...)
gia = FastIsostasy.Simulation(..., model)   # or Gaia.Simulation(...)
sgh = Ypopagos.Simulation(..., model)
bsl = FastIsostasy.PiecewiseConstantBSL()
tpo = Topios.Simulation(..., model)
tools = Kryonomos.Tools(...)

for t in tspan[1]:dt_coupling:tspan[2]
    step!(ice)
    
    couple1way!(smb, ice)
    diagnose!(smb)
    couple1way!(ice, smb)
    # same for ioi and ghf (unless they become prognostic)

    couple1way!(gia, ice)
    step!(gia)
    couple1way!(ice, gia)
    # same for subglacial hydrology, barystatic sea level and landscape evolution
end

function couple1way!(
    gia::FastIsostasy.Simulation,
    ice::Pagos.Simulation,
    tools::Kryonomos.Tools,
)
    conservative_interpolation!(gia.now.H_ice, ice.H_ice, tools)
end

#=
# Pagos

Includes:
- Thermodynamics
- Material
- Dynamics
- Topography
- Boundary conditions
=#

struct Simulation{I1, I2, SM}
    thermodynamic_integrator::I1    # <: DEIntegrator
    topography_integrator::I2       # <: DEIntegrator
    stateful_model::SM              # <: StatefulModel
end

struct StatefulModel{F, S1, S2, M, T, O, P}
    t::F        # <: AbstractFloat
    now::S1     # <: AbstractState
    ref::S2     # <: AbstractState
    model::M    # <: AbstractModel
    tools::T    # <: AbstractTools
    opts::O     # <: Options
    # params::P   # <: AbstractParameters, actually not needed if we do everything right
                # everything will be stored in the model
end

struct State <: AbstractState
    thermodynamics
    material
    dynamics
    topography
    boundary
end
model = Model(
    thermodynamics = Enthalpy(atol = 1f-4),
    material = Arrhenius(),
    dynamics = Stokes(),
)
function update_dynamics!(dyn, topo, mat, dynamics::DIVA)
    genral_stuff()
    update_dynamics!(..., dynamics_solver)
end
function update_dynamics!(dyn, topo, mat, dynamics::DIVA, ..)
end

function (dyn, topo, mat, dynamics::DIVA, dynamics_solver::RobinsonSwier)
    dynamics_solver.nn_inference!(dyn.vel, )
end

struct ThermodynamicsModel
    solver
end

struct Model
    thermodynamics          # Enthalpy, Temperature <: AbstractThermodynamics
    material                # Arrhenius, ... <: AbstractMaterial
    dynamics                # Stokes, BlatterPattyn, DIVA, HybridSXA, L1L2,
                            # SIA, SSA <: AbstractDynamics
    dynamics_solver         # LinearSolve, PT, DynNet <: AbstractDynamicsSolver
    topography              # <: AbstractTopography
    topography_solver       # Classic, LevelSetFunction <: AbstractTopographySolver
    calving                 # EigenCalving, NoCalving, ... <: AbstractCalving
    friction                # CoulombFriction, WeertmanFriction, ... <: AbstractFriction
end

function step!(sim::Simulation, dt)
    # Update the mass balance
    step!(sim.topography_integrator, dt)

    # Update the thermodynamics
    step!(sim.thermodynamic_integrator, dt)

    # Update the material properties
    diagnose!(
        sim.stateful_model.now.material,
        sim.stateful_model.now.thermodynamics,
        sim.stateful_model.model.material)

    # Update the dynamics
    diagnose!(
        sim.stateful_model.now.dynamics,
        sim.stateful_model.now.material,
        sim.stateful_model.model.dynamics)

    # Update the time
    sim.t += dt
end

#=

A key step is to build the integrators in a general way that is compatible
with OrdinaryDiffEq.jl, so that we can use optimized & independent time stepping.

=#

function thermodynamic_rate!(dTdt, T, stateful_model, tspan)
    ...
end

function topography_rate!(dTdt, T, stateful_model, tspan)
    ...
end

thermodynamic_integrator = ODEProblem(
    thermodynamic_rate!,
    stateful_model.now.thermodynamics.T,
    stateful_model,
    stateful_model.t,
)

topography_integrator = ODEProblem(
    topography_rate!,
    stateful_model.now.topography.H_ice,
    stateful_model,
    stateful_model.t,
)
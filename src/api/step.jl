"""
$(TYPEDSIGNATURES)

Advance the simulation by one macro-step `Δt_sync`.

`Δt_sync` is the synchronisation interval shared by all components. Each component's
adaptive integrator takes as many internal sub-steps as needed to reach `t + Δt_sync`.
Components are advanced in operator-splitting order; `couple!` calls propagate updated
fields between components before each dependent step.
"""
function step!(sim::Simulation, Δt_sync)
    step!(sim.ice_sheet.topography, sim.ice_sheet.topography.time_stepper, Δt_sync)

    couple!(sim.ice_sheet.dynamics, sim.ice_sheet.topography)
    couple!(sim.ice_sheet.thermodynamics, sim.ice_sheet.topography)

    step!(sim.ice_sheet.dynamics, sim.ice_sheet.dynamics.time_stepper, Δt_sync)

    couple!(sim.ice_sheet.thermodynamics, sim.ice_sheet.dynamics)
    couple!(sim.ice_sheet.material, sim.ice_sheet.dynamics)
    couple!(sim.ice_sheet.topography, sim.ice_sheet.dynamics)

    step!(
        sim.ice_sheet.thermodynamics,
        sim.ice_sheet.thermodynamics.time_stepper,
        Δt_sync,
    )

    couple!(sim.ice_sheet.material, sim.ice_sheet.thermodynamics)

    step!(sim.ice_sheet.material, sim.ice_sheet.material.time_stepper, Δt_sync)

    couple!(sim.ice_sheet.dynamics, sim.ice_sheet.material)

    step!(sim.timer, Δt_sync)
    is_write_step(sim.timer) && write!(sim.io, sim.ice_sheet)
end

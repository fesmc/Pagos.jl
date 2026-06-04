"""
$(TYPEDSIGNATURES)

Advance the ice sheet model until time `t`.
"""
function step!(icesheet::IceSheet, t)
    # -- Update the velocity field --
    step!(icesheet, icesheet.dynamics, t)

    # -- Update the thickness field --
    step!(icesheet, icesheet.topography, t)

    # -- Update the temperature field --
    step!(icesheet, icesheet.thermodynamics, t)

    # -- Update the calving front --
    calving!(icesheet, icesheet.topography)

    basal_friction!(icesheet, icesheet.dynamics)
    effective_pressure!(icesheet, icesheet.dynamics)

end
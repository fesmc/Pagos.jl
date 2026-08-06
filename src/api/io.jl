struct OutputWriter
    logging::Any     # for progressmeter, timing etc.
    netcdf::Any      # for NetCDF output
    restart::Any     # for restart output (jld2)
end

# @dev TODO: design a lazy output path — variables computed on-the-fly at netcdf write
# time rather than stored in memory. Sketch:
#
#   netcdf = NetCDFOutputWriter(
#       "output.nc",
#       save = [:thickness, :velocity, :driving_stress, :mass_fluxes],
#       diagnose = [:driving_stress, :mass_fluxes],
#   )
#
# `save` lists variables written to the netcdf file; `diagnose` lists ones computed
# on-the-fly at write time. Consider collapsing both into one list with a per-variable
# flag for how it's computed, instead of two separate lists.
abstract type AbstractOutput end
struct OutputWriter
    logging::Any     # for progressmeter, timing etc.
    netcdf::Any      # for NetCDF output
    restart::Any     # for restart output (jld2)
end

"""
I would like to define a list of variables that are only computed when writing out the netcdf file. They should not be stored in memory, but computed on-the-fly. In terms of API, this should look like this:

netcdf = NetCDFOutputWriter(
    "output.nc",
    save = [
        :thickness,
        :velocity,
        :driving_stress,
        :mass_fluxes,
    ],
    diagnose = [
        :driving_stress,
        :mass_fluxes,
    ],
)

`save` provides the list of variables that are written to the netcdf file, while `diagnose` provides the list of variables that are computed on-the-fly when writing out the netcdf file. 

I wonder whether we could collapse both very easily into a single list of variables, and then have a flag in our internal dictionary that indicates how the variable is computed (on-the-fly or stored in memory).
"""
abstract type AbstractOutput end
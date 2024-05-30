using NCDatasets
include("PTloop.jl")

Lx = 6000e3
power_vec, dx, nx = sample_resolutions((1, 5), Lx)
dtau = nx ./ nx[2]
slab = TestSlab()
runtime = rand(length(dx))

for i in eachindex(dx)
    icesheet = slab_icesheet(slab, dx[i]; nx = nx[i], ny = nx[i], dtau_scaling = dtau[i])
    pseudo_transient!(icesheet)
    @show extrema(icesheet.state.ux)
    pseudodotvel_time, pseudovel_time, solve_time = time_slab_problem(icesheet)
    @show solve_time
    println("------------------")
    runtime[i] = solve_time
end

NCDataset("PTloop-runtimes.nc", "c") do ds
    defDim(ds, "resolution", length(runtime))
    defVar(ds, "resolution", dx, ("resolution",), attrib = Dict(
        "units" => "m",
        "comments" => "Resolution of slab experiment",
    ))
    defVar(ds, "runtime", runtime, ("resolution",), attrib = Dict(
           "units" => "seconds",
           "comments" => "Minimum runtime of pseudo-transient loop for given resolution",
    ))
end
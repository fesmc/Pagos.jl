include("PTloop.jl")
Lx, dx = 6000e3, 16e3
nx = round(Int, Lx / dx)
slab = TestSlab()
icesheet = slab_icesheet(slab, dx; nx = nx, ny = nx, dtau_scaling = 1.0)
pseudo_transient!(icesheet)
pseudodotvel_time, pseudovel_time, solve_time = btime_slab_problem(icesheet)
# (0.005976416, 4.84e-5, 0.008410827)

# In comparison, the linear solver gives: 3.573 s (1066 allocations: 1010.89 MiB)
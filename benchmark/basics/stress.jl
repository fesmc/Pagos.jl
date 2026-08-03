using Pagos
using Printf
import KernelAbstractions

# RegularGrid fixes nz = 1 (2D). Rebuild with nz > 1 for the full-column (3D) case.
function grid3d(T, lx, ly, dx, dy, nz)
    g = RegularGrid(T, lx, ly, dx, dy)
    z = collect(range(T(0), T(1); length = nz))
    return RegularGrid(g.nx, g.ny, nz, g.x, g.y, z, g.dx, g.dy, g.dz,
        g.Lon, g.Lat, g.area, g.distortion, g.basins, g.regions)
end

function fill_inputs!(mech, mat)
    fill!(mat.eta_ice, 2.0)
    fill!(mech.strain_rate_dxx,  3.0); fill!(mech.strain_rate_dyy, -1.0)
    fill!(mech.strain_rate_dzz, -2.0); fill!(mech.strain_rate_dxy,  0.5)
    fill!(mech.strain_rate_dxz,  1.0); fill!(mech.strain_rate_dyz, -0.5)
    return nothing
end

# min-of-N timer (kernels are ms-scale, well above @elapsed resolution).
# Operators launch asynchronously on GPU, so the timed region must synchronize —
# otherwise @elapsed measures launch overhead, not kernel execution.
function bench(f, args...; backend = KernelAbstractions.CPU(), samples = 300)
    f(args...)                       # warmup / compile
    KernelAbstractions.synchronize(backend)
    best = Inf
    for _ in 1:samples
        best = min(best, @elapsed begin
            f(args...)
            KernelAbstractions.synchronize(backend)
        end)
    end
    return best
end

function run_case(name, grid)
    mech = MechanicState(grid)
    mat  = MaterialState(grid)
    fill_inputs!(mech, mat)

    backend = KernelAbstractions.get_backend(mat.eta_ice)
    t = bench(deviatoric_stress!, mech, mat; backend)
    n = length(mech.stress_xx)
    @printf("%-20s N=%9d   deviatoric_stress! = %7.3f ms\n", name, n, 1e3t)
    return nothing
end

println("Julia threads = ", Threads.nthreads())
println("fused KernelAbstractions kernel (single pass)\n")

run_case("2D 512x512",    RegularGrid(Float64, 511.0, 511.0, 1.0, 1.0))
run_case("2D 1024x1024",  RegularGrid(Float64, 1023.0, 1023.0, 1.0, 1.0))
run_case("2D 2048x2048",  RegularGrid(Float64, 2047.0, 2047.0, 1.0, 1.0))
run_case("3D 256x256x20",  grid3d(Float64, 255.0, 255.0, 1.0, 1.0, 20))
run_case("3D 512x512x20",  grid3d(Float64, 511.0, 511.0, 1.0, 1.0, 20))

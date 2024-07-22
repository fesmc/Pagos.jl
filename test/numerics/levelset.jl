using LinearAlgebra
using CairoMakie

# Initialize the level set function
function initialize_level_set(grid, center, radius)
    X, Y, Z = grid
    return sqrt.((X .- center[1]).^2 .+ (Y .- center[2]).^2 .+ (Z .- center[3]).^2) .- radius
end

# Define the velocity field
function define_velocity_field(grid)
    # Example: uniform velocity field
    u = ones(size(grid[1]))
    v = zeros(size(grid[1]))
    w = zeros(size(grid[1]))
    return u, v, w
end

# Upwind derivative
function upwind_derivative(phi, axis, direction)
    if direction == "forward"
        return circshift(phi, -1, axis) .- phi
    else
        return phi .- circshift(phi, 1, axis)
    end
end

# Update the level set function
function update_level_set(phi, u, v, w, dt, dx)
    phi_x_f = upwind_derivative(phi, 1, "forward")
    phi_x_b = upwind_derivative(phi, 1, "backward")
    phi_y_f = upwind_derivative(phi, 2, "forward")
    phi_y_b = upwind_derivative(phi, 2, "backward")
    phi_z_f = upwind_derivative(phi, 3, "forward")
    phi_z_b = upwind_derivative(phi, 3, "backward")

    phi_t = min.(u, 0) .* phi_x_f .+ max.(u, 0) .* phi_x_b .+
            min.(v, 0) .* phi_y_f .+ max.(v, 0) .* phi_y_b .+
            min.(w, 0) .* phi_z_f .+ max.(w, 0) .* phi_z_b

    return phi .- dt * phi_t / dx
end

# Reinitialize the level set function
function reinitialize_level_set(phi, dx)
    # Placeholder for reinitialization step
    return phi
end

# Grid setup
Nx, Ny, Nz = 100, 100, 100
x = range(-1, 1, length=Nx)
y = range(-1, 1, length=Ny)
z = range(-1, 1, length=Nz)
X, Y, Z = [reshape(v, Nx, 1, 1) for v in x], [reshape(v, 1, Ny, 1) for v in y], [reshape(v, 1, 1, Nz) for v in z]

# Initialization
phi = initialize_level_set((X, Y, Z), [0, 0, 0], 0.5)
u, v, w = define_velocity_field((X, Y, Z))

# Time-stepping parameters
dt = 0.01
dx = x[2] - x[1]
num_time_steps = 100
reinit_interval = 10

# Time-stepping loop
for n in 1:num_time_steps
    phi = update_level_set(phi, u, v, w, dt, dx)
    if n % reinit_interval == 0
        phi = reinitialize_level_set(phi, dx)
    end
end

# Result visualization
# contourf(x, y, phi[:, :, Nz ÷ 2], levels = [-0.5, 0, 0.5])
# title!("Level Set Function at Final Time Step")
# xlabel!("x")
# ylabel!("y")

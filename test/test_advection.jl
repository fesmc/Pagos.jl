# Here is a simple test script to see how to make an explicit advection solver using ODEProblem

cd(@__DIR__)
import Pkg; Pkg.activate("../")
using Revise

using DifferentialEquations
using CairoMakie

struct MassConservationParams
    ux::Matrix
    uy::Matrix
    a::Matrix
    nx::Integer
    ny::Integer
    dx::Float64
    dy::Float64
end

# RHS function
function rhs!(dH_vec, H_vec, p::MassConservationParams, t)

    # Get variables in 2D form
    H = reshape(H_vec, p.nx, p.ny)
    dH = zeros(p.nx, p.ny)

    for i in 2:p.nx-1, j in 2:p.ny-1
        # Upwind fluxes (ux > 0, uy > 0)
        flux_x = p.ux[i,j] * (H[i,j] - H[i-1,j]) / p.dx
        flux_y = p.uy[i,j] * (H[i,j] - H[i,j-1]) / p.dy
        dH[i,j] = - (flux_x + flux_y) + p.a[i,j]
    end

    # Neumann boundary conditions (copy nearest interior values)
    dH[1,:] .= dH[2,:]
    dH[end,:] .= dH[end-1,:]
    dH[:,1] .= dH[:,2]
    dH[:,end] .= dH[:,end-1]

    dH_vec[:] .= vec(dH)
end

# Visualization at time t
function plot_H(t)
    H = reshape(sol(t), nx, ny)
    fig, ax, hm = heatmap(x, y, H,colorrange=(0,1))
    ax.title ="Ice Thickness at t = $t"
    fig
end

# Grid parameters
nx, ny = 30, 50
Lx, Ly = 1.0, 1.0
x = LinRange(0, Lx, nx)
y = LinRange(0, Ly, ny)
dx = x[2] - x[1]
dy = y[2] - y[1]

# Velocity field: constant rightward and upward flow
ux = fill(0.2, nx, ny)
uy = fill(0.1, nx, ny)

# Source term a(x, y)
a = [0.1 * sin(2π * x[i]) * sin(2π * y[j]) for i in 1:nx, j in 1:ny]

# Initial condition: Gaussian bump in center
H0 = [exp(-100 * ((x[i] - 0.5)^2 + (y[j] - 0.5)^2)) for i in 1:nx, j in 1:ny]
H0_vec = vec(H0)  # flatten into 1D vector

# Define input variables for mass conservation problem
params = MassConservationParams(ux, uy, a, nx, ny, dx, dy)

# Problem setup
tspan = (0.0, 10.0)
prob = ODEProblem(rhs!, H0_vec, tspan, params)

sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)

# Plot at selected time points
plot_H(0.0)
plot_H(0.25)
plot_H(0.5)
plot_H(1.0)
plot_H(5.0)
plot_H(10.0)
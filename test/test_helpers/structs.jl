include("utils.jl")

"""
    TestSlab

A class representing a test slab.
"""
struct TestSlab
    H
    mu
    beta
    alpha
    rho
    g
    eta
    F2
    ub
    u
end

function TestSlab(; H = 1000.0, mu = 1e5, beta = 1e3, alpha = 1e-3, rho = 910.0, g = 9.81)
    eta = f_eta(beta, H, mu)
    F2 = f_F2(H, mu)
    ub = f_ub(rho, g, H, alpha, beta)
    u = f_u(rho, g, H, alpha, beta, mu)
    return TestSlab(H, mu, beta, alpha, rho, g, eta, F2, ub, u)
end

"""
    TestStream

A class representing a test stream following Schoof (2006).
"""
struct TestStream
    H
    alpha
    W
    L
    m
    rho
    g
    n_glen
    rf
    yc
    ux
    tau_c
    c_bed
    beta
end

function TestStream(; dx = 1e3, H = 1e3, rf = 1e-16, alpha = 1e-3, rho = 910.0, g = 9.81, n_glen = 3.0,
    W = 25e3, m = 1.55, xmax = 140e3)
    
    T = Float64
    lx = xmax
    ly = 4 * W
    dy = dx
    domain = Domain(T, lx, ly, dx, dy)
end
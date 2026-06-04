module PagosCUDSSExt

using Pagos, CUDA, CUDSS

# Direct CUDSS solve for GPU-backed LinearDynamicsSolver2D (CSR matrix required).
# Bypasses LinearSolve entirely, so the LinearSolveCUDSSExt version conflict
# is irrelevant.
#
# The sparsity pattern is fixed for a given grid, so symbolic analysis only
# runs once. Numeric factorization runs on every call (matrix values change).
function Pagos.velocity!(lsd::Pagos.LinearDynamicsSolver2D{<:Any,<:Any,<:CuVector}; kw...)
    solver = if lsd.solver_cache[] === nothing
        s = CudssSolver(lsd.A, "G", 'F')
        cudss("analysis", s, lsd.u, lsd.b)
        lsd.solver_cache[] = s
        s
    else
        lsd.solver_cache[]::CudssSolver
    end
    cudss("factorization", solver, lsd.u, lsd.b)
    cudss("solve",         solver, lsd.u, lsd.b)
    return nothing
end

end

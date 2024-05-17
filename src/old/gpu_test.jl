using SparseArrays, LinearAlgebra
using CUDA, LinearSolve, Random
CUDA.allowscalar(false)

n = 400
Is = collect(1:n^2)[rand(n^2) .> 0.95]
Js = shuffle(Is)
Vs = rand(length(Is))
A = sparse(Is, Js, Vs) + I
b = rand(size(A, 1))

A_gpu = CUDA.CUSPARSE.CuSparseMatrixCSR(A)
b_gpu = CuVector(b)
x_gpu = CUDA.zeros(Float64, length(b))

@btime CUDA.CUSOLVER.csrlsvqr!($A_gpu, $b_gpu, $x_gpu, 1e-4, one(Cint), 'O')
@btime x = $A \ $b

prob = LinearProblem(A, b)
@btime sol = solve(prob)
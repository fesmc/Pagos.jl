A = [1.0 2.0; 3.0 4.0]
println(hasnan(A)) # false

A = [1.0 2.0; 3.0 NaN]
println(hasnan(A)) # true
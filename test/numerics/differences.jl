n = 191 * 2
Z = (1:n) * ones(n)'
Z_dx = similar(Z)
nx, ny = size(Z)
@btime dx!($Z_dx, $Z, 1.0, nx)
@btime dx($Z, 1.0)
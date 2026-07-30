"""
    hasnan(x::AbstractArray)

Check if an array has NaN values. Accepts a Chmy `Field` too, and then inspects its
interior only: the reduction goes through [`asarray`](@ref), so it runs on the underlying
array on every backend rather than falling back to a `Field`'s scalar `getindex`. Halo
cells are excluded, which is what you want — an unfilled halo is not a NaN in the solution.
"""
hasnan(x::AbstractArray) = any(isnan, asarray(x))
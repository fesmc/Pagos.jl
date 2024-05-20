"""
    hasnan(x::AbstractArray)

Check if an array has NaN values.
"""
hasnan(x::AbstractArray) = any(isnan, x)
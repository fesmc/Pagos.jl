"""
$(TYPEDSIGNATURES)

Saturate a value `x` between `xmin` and `xmax`.
"""
function saturate(x, xmin, xmax)
    return min(max(x, xmin), xmax)
end

## Using masks

Based on masks, `CartesianIndices` can be created and allow iterating only over a
desired region. For instance:

```julia
icemask = H .> 0
indices2D = CartesianIndices((nx, ny))
icemask_indices = indices2D[icemask]
pseudo_transient!(args..., icemask_indices)
```

will be faster than its equivalent without. The `icemask_indices` should be created
only once at each update of the ice thickness.
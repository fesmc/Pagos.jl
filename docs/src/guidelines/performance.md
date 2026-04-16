# Performance

## Systematical benchmark

In future, the documentation will include some benchmarks. Every PR that modifies existing routines needs to prove that the code is not slowed down or, if so, for a good reason.

## Using `map!`

The in-place `map!` applies functions to arrays without allocation. It combines well with different array types, e.g., `SubArray`, `CuArray`, etc. It might be less performant for specific cases but is much more general and should therefore be preferred. Simple benchmarks showed that applying `map!` over the full array is faster than applying it only to a view of the array. Therefore, if some operation is not defined outside of a mask, it is better to apply `map!` to the full array and use a conditional inside the function.
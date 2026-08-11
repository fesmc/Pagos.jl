"""
$(TYPEDSIGNATURES)

Generic elementwise-map kernel shared by every pointwise physics function
(`creep!`, `calving_rate!`, `effective_pressure!`, ...): writes
`out[I] = f(args[I]..., extra...)` at every index `I`, where `args` are the
per-cell array arguments and `extra` are non-indexed trailing arguments
(parameter structs, scalars such as grid spacing).

!!! warning "Do not reassign the offset-adjusted index before closing over it"
    `Iglob` must be bound once and not reassigned before the `map` closure
    below captures it. Reassigning a captured variable (e.g. `I = I + O`, then
    closing over `I`) boxes it (`Core.Box`), which is silently slower on CPU
    but fails GPU compilation outright ("unsupported dynamic function
    invocation") — see `pagos-roadmap/chmy.md`, appendix, 2026-07-31.
"""
@kernel inbounds = true function _pointwise!(f::F, out, args, extra, O) where {F}
    Ilocal = @index(Global, NTuple)
    Iglob = Ilocal + O
    out[Iglob...] = f(map(a -> a[Iglob...], args)..., extra...)
end

"""
$(TYPEDSIGNATURES)

Launch [`_pointwise!`](@ref) over the plain array `out`: `out[I] = f(args[I]...,
extra...)` at every index `I` of `out`. `args` and `extra` are always tuples,
even for a single array/parameter argument. Backend is read off `out`
(`KernelAbstractions.get_backend`); no grid or `Runtime` is involved.
"""
function pointwise!(f::F, out::AbstractArray, args::Tuple, extra::Tuple = ()) where {F}
    backend = KernelAbstractions.get_backend(out)
    _pointwise!(backend)(f, out, args, extra, Offset(); ndrange = size(out))
    return nothing
end

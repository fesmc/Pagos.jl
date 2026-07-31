using Pagos
using JET
using Test

@testset "JET static analysis" begin
    # `:typo` mode reports undefined bindings, invalid field accesses and unmatched
    # keyword arguments. It deliberately skips JET's `:basic`/`:sound` dynamic-dispatch
    # and runtime-dispatch diagnostics, which are extremely noisy on Tullio- and
    # KernelAbstractions-generated code and are not actionable here.
    #
    # `target_modules = (Pagos,)` scopes the *reported* frames to Pagos's own source,
    # so pre-existing false positives inside Base/LinearAlgebra/KernelAbstractions
    # internals (triggered merely by calling generic dense/sparse linear algebra or
    # GPU-kernel machinery) don't show up here.
    result = JET.report_package(
        Pagos;
        mode = :typo,
        target_modules = (Pagos,),
        toplevel_logger = nothing,
    )
    reports = JET.get_reports(result)

    # `TillEffectivePressure()` (called with no type parameter) can't infer `T` for its
    # `overburden::OverburdenEffectivePressure{T}` field default. Nothing in Pagos calls
    # this bare form — every other construction of `Constants`/`*EffectivePressure`
    # explicitly specifies the type parameter, e.g. `TillEffectivePressure{Float64}()` —
    # so it is tracked here instead of fixed outright (the recursive type-parameter
    # default would need a small redesign of the struct).
    # TODO: remove once `TillEffectivePressure()` gets a real default `T`, or the bare
    # no-type-param constructor is dropped from the public API.
    is_known_issue(r) = r isa JET.UndefVarErrorReport &&
        string(r.var) == "Pagos.T" &&
        any(f -> occursin("effective_pressure.jl", String(f.file)), r.vst)

    unexpected = filter(!is_known_issue, reports)
    if !isempty(unexpected)
        show(result)
        println()
    end
    @test isempty(unexpected)
end

using Pagos
using Aqua
using Test

@testset "Aqua.jl quality assurance" begin
    Aqua.test_all(
        Pagos;
        # `undefined_exports` and `undocumented_names` are relaxed until the dead exports
        # (names moved to `src/legacy/` or still WIP, e.g. `ζ_aa`, `AbstractFloatingCalving`,
        # `depthaveraged_velocity!`, `scaledstrainrate!`)
        # are either removed from the export lists in `src/Pagos.jl` or given definitions.
        # Once that is done, both can be re-enabled: every *defined* public name is already
        # documented. See also the `undocumented_names`/`undefined_exports` findings.
        undefined_exports = false,
        undocumented_names = false,
        # TODO: `FastGaussQuadrature` and `FillArrays` are currently only used from
        # un-`include`d `src/legacy/` code, and `NCDatasets` only from the test suite / I/O
        # helpers. Remove the ignore list once they are wired into `src/` (or dropped as deps).
        stale_deps = (ignore = [:FastGaussQuadrature, :NCDatasets, :FillArrays],),
    )
end

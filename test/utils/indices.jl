using Pagos
using Test

@testset "AbstractIndexing" begin
    i1, i2 = 1, 10
    strict     = StrictIndexing(i1, i2)
    flat       = FlatIndexing(i1, i2)
    reflective = ReflectiveIndexing(i1, i2)
    periodic   = PeriodicIndexing(i1, i2)

    @test index(2, -1, strict)     == 1
    @test index(1, -3, flat)       == 1
    @test index(10, 1, reflective) == 9
    @test index(1, -1, periodic)   == 10
end

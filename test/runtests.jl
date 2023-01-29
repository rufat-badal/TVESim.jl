using TVESim
using Test
using Random
using JuMP

const RNG_SEED = 42

@testset "Linear algebra operations" begin
    rng = MersenneTwister(RNG_SEED)
    model = Model()

    m, n = 100, 120
    scalar = 23.0
    A_val = rand(rng, m, n)
    @NLexpression(model, A[i=1:m, j=1:n], A_val[i, j])

    @test scalar * A_val ≈ value.(scalar * A)
end

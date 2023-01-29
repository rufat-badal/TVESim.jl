using TVESim
using Test
using Random
using JuMP

const RNG_SEED = 42

@testset "Linear algebra operations" begin
    rng = MersenneTwister(RNG_SEED)
    n = 100
    model = Model()

    A_val = rand(rng, n, n)
    @NLexpression(model, A[i=1:n, j=1:n], A_val[i, j])

    @test 2A_val ≈ value.(2A)
end

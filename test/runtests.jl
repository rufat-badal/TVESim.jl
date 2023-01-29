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
    A = @NLexpression(model, [i = 1:m, j = 1:n], A_val[i, j])

    @test value.(scalar * A) ≈ scalar * A_val
    @test value.(-A) ≈ -A_val

    B_val = rand(rng, m, n)
    B = @NLexpression(model, [i = 1:m, j = 1:n], B_val[i, j])

    @test value.(A + B) ≈ A_val + B_val
end

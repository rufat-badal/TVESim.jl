using TVESim
using Test
using Random
using JuMP
using LinearAlgebra

const RNG_SEED = 42

function get_matrix_pair(rng, model, m, n)
    A_val = rand(rng, m, n)
    A = @NLexpression(model, [i = 1:m, j = 1:n], A_val[i, j])
    A, A_val
end

@testset "Linear algebra operations" begin
    rng = MersenneTwister(RNG_SEED)
    model = Model()

    m, n = 100, 120
    scalar = 23.0
    A, A_val = get_matrix_pair(rng, model, m, n)
    @test value.(scalar * A) ≈ scalar * A_val
    @test value.(-A) ≈ -A_val

    B, B_val = get_matrix_pair(rng, model, m, n)
    @test value.(A + B) ≈ A_val + B_val
    @test value.(A - B) ≈ A_val - B_val

    m, n, l = 100, 123, 78
    A, A_val = get_matrix_pair(rng, model, m, n)
    B, B_val = get_matrix_pair(rng, model, n, l)
    @test value.(A * B) ≈ A_val * B_val 

    n = 100
    A, A_val = get_matrix_pair(rng, model, n, n)
    @test value(tr(A)) ≈ tr(A_val)
end

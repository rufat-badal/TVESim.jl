using TVESim
using Test
using Random
using JuMP
using LinearAlgebra

const RNG_SEED = 42

function get_matrix_pair(rng, model, m, n)
    A_val = rand(rng, m, n)
    A = [TVESim.AdvancedNonlinearExpression(model, @NLexpression(model, A_val[i, j])) for i in 1:m, j in 1:n]
    A, A_val
end

@testset "Linear algebra operations" begin
    rng = MersenneTwister(RNG_SEED)
    model = Model()

    m, n = 100, 120
    λ_val = 23
    λ = TVESim.AdvancedNonlinearExpression(model, λ_val)
    A, A_val = get_matrix_pair(rng, model, m, n)
    @test value.(λ_val * A) == λ_val * A_val
    @test value.(λ * A) == λ_val * A_val
    @test value.(λ.expression * A) == λ_val * A_val
    @test value.(A * λ.expression) == λ_val * A_val
    @test value.(A * λ) == λ_val * A_val
    @test value.(A / λ_val) == A_val / λ_val
    @test value.(A / λ) == A_val / λ_val
    @test value.(A / λ.expression) == A_val / λ_val
    @test value.(-A) == -A_val
    @test value.(TVESim.transpose(A)) == transpose(A_val)
    @test value(TVESim.norm_sqr(A)) ≈ sum(A_val .^ 2)

    B, B_val = get_matrix_pair(rng, model, m, n)
    @test value.(A + B) == A_val + B_val
    @test value.(A_val + B) == A_val + B_val
    @test value.(A + B_val) == A_val + B_val
    @test value.(A - B) == A_val - B_val
    @test value.(A_val - B) == A_val - B_val
    @test value.(A - B_val) == A_val - B_val
    @test value(TVESim.dot(A, B)) ≈ dot(A_val, B_val)
    @test value(TVESim.dot(A_val, B)) ≈ dot(A_val, B_val)
    @test value(TVESim.dot(A, B_val)) ≈ dot(A_val, B_val)

    m, n, l = 100, 123, 78
    A, A_val = get_matrix_pair(rng, model, m, n)
    B, B_val = get_matrix_pair(rng, model, n, l)
    @test value.(A * B) ≈ A_val * B_val
    @test value.(A_val * B) ≈ A_val * B_val
    @test value.(A * B_val) ≈ A_val * B_val

    n = 100
    A, A_val = get_matrix_pair(rng, model, n, n)
    @test value(TVESim.tr(A)) ≈ tr(A_val)

    n = 7
    A, A_val = get_matrix_pair(rng, model, n, n)
    @test value(TVESim.det(A)) ≈ det(A_val)
    @test value.(A * TVESim.adjugate(A)) ≈ det(A_val) * Matrix(I, n, n)
end

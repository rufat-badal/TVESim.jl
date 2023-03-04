using TVESim
using Test
using Random
using JuMP
using LinearAlgebra

const RNG_SEED = 42

function get_matrix_pair(rng, model, m, n)
    X_val = rand(rng, m, n)
    X = @NLparameter(model, [i = 1:m, j = 1:n] == X_val[i, j])
    TVESim.jumpexpression_array(model, X), X_val
end

# rng = MersenneTwister(RNG_SEED)
# model = Model()
# m, n = 2, 2
# X, X_val = get_matrix_pair(rng, model, m, n)
# λ_val = 3
# λ = TVESim.JuMPExpression(model, @NLparameter(model, value = λ_val))

@testset "Linear algebra operations" begin
    rng = MersenneTwister(RNG_SEED)
    model = Model()

    m, n = 100, 120
    λ_val = 23
    λ = TVESim.JuMPExpression(model, @NLparameter(model, value = λ_val))
    X, X_val = get_matrix_pair(rng, model, m, n)
    @test value.(X) == X_val
    @test value.(-X) == -X_val
    @test value.(λ_val * X) == λ_val * X_val
    @test value.(X * λ_val) == X_val * λ_val
    @test value.(λ * X) == λ_val * X_val
    @test value.(X * λ) == X_val * λ_val
    @test value.(λ.expr * X) == λ_val * X_val
    @test value.(X * λ.expr) == X_val * λ_val
    @test value.(X / λ_val) == X_val / λ_val
    @test value.(X / λ) == X_val / λ_val
    @test value.(X / λ.expr) == X_val / λ_val
    @test value.(transpose(X)) == transpose(X_val)
    @test value(sum(X .^ 2)) ≈ sum(X_val .^ 2)

    Y, Y_val = get_matrix_pair(rng, model, m, n)
    @test value.(X + Y) == X_val + Y_val
    @test value.(X_val + Y) == X_val + Y_val
    @test value.(X + Y_val) == X_val + Y_val
    @test value.(X - Y) == X_val - Y_val
    @test value.(X_val - Y) == X_val - Y_val
    @test value.(X - Y_val) == X_val - Y_val
    @test value(dot(X, Y)) ≈ dot(X_val, Y_val)
    @test value(dot(X_val, Y)) ≈ dot(X_val, Y_val)
    @test value(dot(X, Y_val)) ≈ dot(X_val, Y_val)

    m, n, l = 100, 123, 78
    X, X_val = get_matrix_pair(rng, model, m, n)
    Y, Y_val = get_matrix_pair(rng, model, n, l)
    @test value.(X * Y) ≈ X_val * Y_val
    @test value.(X_val * Y) ≈ X_val * Y_val
    @test value.(X * Y_val) ≈ X_val * Y_val

    n = 100
    X, X_val = get_matrix_pair(rng, model, n, n)
    @test value(tr(X)) ≈ tr(X_val)

    n = 7
    X, X_val = get_matrix_pair(rng, model, n, n)
    @test value(TVESim.det(X)) ≈ det(X_val)
    @test value.(X * TVESim.adjugate(X)) ≈ det(X_val) * Matrix(I, n, n)
end

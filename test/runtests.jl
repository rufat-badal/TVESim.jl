using TVESim
using Test
using Random
using JuMP
using LinearAlgebra

const RNG_SEED = 42

function JuMP.value(e::TVESim.AdvancedNonlinearExpression)
    value(e._expression)
end

function get_matrix_pair(rng, model, m, n)
    A_val = rand(rng, m, n)
    A = [TVESim.AdvancedNonlinearExpression(model, @NLexpression(model, A_val[i, j])) for i in 1:m, j in 1:n]
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
    @test value.(transpose(A)) ≈ transpose(A_val)
    @test value(sum(A .^ 2)) ≈ sum(A_val .^ 2)

    B, B_val = get_matrix_pair(rng, model, m, n)
    @test value.(A + B) ≈ A_val + B_val
    # @test value.(A - B) ≈ A_val - B_val
    # @test value.(dot(A, B)) ≈ dot(A_val, B_val)

#     # m, n, l = 100, 123, 78
#     # A, A_val = get_matrix_pair(rng, model, m, n)
#     # B, B_val = get_matrix_pair(rng, model, n, l)
#     # @test value(A * B) ≈ A_val * B_val

#     # n = 100
#     # A, A_val = get_matrix_pair(rng, model, n, n)
#     # @test value(tr(A)) ≈ tr(A_val)

#     # n = 7
#     # A, A_val = get_matrix_pair(rng, model, n, n)
#     # @test value(det(A)) ≈ det(A_val)
#     # @test value(A * TVESim.adjugate(A)) ≈ det(A_val) * Matrix(I, n, n)
end

# model = Model()
# e1 = TVESim.AdvancedNonlinearExpression(model, @NLexpression(model, 1.0))
# e2 = TVESim.AdvancedNonlinearExpression(model, @NLexpression(model, 2.0))
# e3 = 3.0
# display(e1)
# display(e2)
# display(-e1)
# display(e1 + e2)
# display(e1 + e3)
# display(e3 + e1)
# display(e1 - e2)
# display(e1 - e3)
# display(e3 - e1)
# display(e1 * e2)
# display(e1 * e3)
# display(e3 * e1)
# value(e3 * e2)

# rng = MersenneTwister(RNG_SEED)
# model = Model()
# m, n = 100, 120
# A, A_val = get_matrix_pair(rng, model, m, n)
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

# @testset "Linear algebra operations" begin
#     rng = MersenneTwister(RNG_SEED)
#     model = Model()

#     m, n = 10, 12
#     scalar = 23.0
#     A, A_val = get_matrix_pair(rng, model, m, n)
#     @test value.(scalar * A) == scalar * A_val
#     @test value.(-A) == -A_val
#     @test value.(transpose(A)) == transpose(A_val)
#     @test value(TVESim.norm_sqr(A)) ≈ sum(A_val .^ 2)

#     B, B_val = get_matrix_pair(rng, model, m, n)
#     @test value.(A + B) == A_val + B_val
#     @test value.(A_val + B) == A_val + B_val
#     @test value.(A + B_val) == A_val + B_val
#     @test value.(A - B) == A_val - B_val
#     @test value.(A_val - B) == A_val - B_val
#     @test value.(A - B_val) == A_val - B_val
#     @test value(dot(A, B)) == dot(A_val, B_val)
#     @test value(dot(A_val, B)) == dot(A_val, B_val)
#     @test value(dot(A, B_val)) == dot(A_val, B_val)

#     m, n, l = 100, 123, 78
#     A, A_val = get_matrix_pair(rng, model, m, n)
#     B, B_val = get_matrix_pair(rng, model, n, l)
#     @test value.(A * B) ≈ A_val * B_val
#     @test value.(A_val * B) ≈ A_val * B_val
#     @test value.(A * B_val) ≈ A_val * B_val

#     n = 100
#     A, A_val = get_matrix_pair(rng, model, n, n)
#     @test value(tr(A)) ≈ tr(A_val)

#     n = 7
#     A, A_val = get_matrix_pair(rng, model, n, n)
#     # @test value(det(A)) ≈ det(A_val)
#     # @test value(A * TVESim.adjugate(A)) ≈ det(A_val) * Matrix(I, n, n)
# end

rng = MersenneTwister(RNG_SEED)
model = Model()
m, n = 4, 5
A, A_val = get_matrix_pair(rng, model, m, n)
B, B_val = get_matrix_pair(rng, model, m, n)
C, C_val = get_matrix_pair(rng, model, m, m)
display(value.(A))
display(value.(TVESim.minor(A, 1, 2)))
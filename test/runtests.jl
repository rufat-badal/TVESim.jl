using TVESim
using Test
using Random
using JuMP
using LinearAlgebra

const RNG_SEED = 42

function JuMP.value(A::TVESim.NLExprMatrix)
    value.(A._matrix)
end

function get_matrix_pair(rng, model, m, n)
    A_val = rand(rng, m, n)
    A = @NLexpression(model, [i = 1:m, j = 1:n], A_val[i, j])
    TVESim.NLExprMatrix(model, A), A_val
end

@testset "Linear algebra operations" begin
    rng = MersenneTwister(RNG_SEED)
    model = Model()

    m, n = 100, 120
    scalar = 23.0
    A, A_val = get_matrix_pair(rng, model, m, n)
    @test value(scalar * A) ≈ scalar * A_val
    @test value(-A) ≈ -A_val
    @test value(transpose(A)) ≈ transpose(A_val)
    # test getindex for NLExprMatrix
    A_internal_matrix = Matrix{NonlinearExpression}(undef, m, n)
    for j in 1:n
        for i in 1:m
            A_internal_matrix[i, j] = A._matrix[i, j]
        end
    end
    @test A_internal_matrix == A._matrix
    @test value(TVESim.norm_sqr(A)) ≈ sum(A_val .^ 2)
    A_tilde = TVESim.NLExprMatrix([A[i, j] for i in 1:m, j in 1:n])
    @test value(A_tilde) == value(A)

    B, B_val = get_matrix_pair(rng, model, m, n)
    @test value(A + B) ≈ A_val + B_val
    @test value(A - B) ≈ A_val - B_val
    @test value(dot(A, B)) ≈ dot(A_val, B_val)

    m, n, l = 100, 123, 78
    A, A_val = get_matrix_pair(rng, model, m, n)
    B, B_val = get_matrix_pair(rng, model, n, l)
    @test value(A * B) ≈ A_val * B_val

    n = 100
    A, A_val = get_matrix_pair(rng, model, n, n)
    @test value(tr(A)) ≈ tr(A_val)

    n = 7
    A, A_val = get_matrix_pair(rng, model, n, n)
    @test value(det(A)) ≈ det(A_val)
    @test value(A * TVESim.adjugate(A)) ≈ det(A_val) * Matrix(I, n, n)
end

# rng = MersenneTwister(RNG_SEED)
# model = Model()
# m, n = 100, 120
# A, A_val = get_matrix_pair(rng, model, m, n)
# B, B_val = get_matrix_pair(rng, model, m, n)
function Base.:*(scalar::Number, A::Matrix{JuMP.NonlinearExpression})
    m, n = size(A)
    model = A[1, 1].model
    A_scaled = Matrix{JuMP.NonlinearExpression}(undef, m, n)
    for j in 1:n
        for i in 1:m
            A_scaled[i, j] = JuMP.@NLexpression(model, scalar * A[i, j])
        end
    end
    A_scaled
end

function Base.:-(A::Matrix{JuMP.NonlinearExpression})
    m, n = size(A)
    model = A[1, 1].model
    minus_A = Matrix{JuMP.NonlinearExpression}(undef, m, n)
    for j in 1:n
        for i in 1:m
            minus_A[i, j] = JuMP.@NLexpression(model, -A[i, j])
        end
    end
    minus_A
end

function checksizematch(A, B)
    size(A) == size(B) || throw(DimensionMismatch("matrix sizes do not match: dimensions are $(size(A)), $(size(B))"))
    return size(A)
end

function Base.:+(A::Matrix{JuMP.NonlinearExpression}, B::Matrix{JuMP.NonlinearExpression})
    m, n = checksizematch(A, B)
    model = A[1, 1].model
    A_plus_B = Matrix{JuMP.NonlinearExpression}(undef, m, n)
    for j in 1:n
        for i in 1:m
            A_plus_B[i, j] = JuMP.@NLexpression(model, A[i, j] + B[i, j])
        end
    end
    A_plus_B
end

function Base.:-(A::Matrix{JuMP.NonlinearExpression}, B::Matrix{JuMP.NonlinearExpression})
    m, n = checksizematch(A, B)
    model = A[1, 1].model
    A_minus_B = Matrix{JuMP.NonlinearExpression}(undef, m, n)
    for j in 1:n
        for i in 1:m
            A_minus_B[i, j] = JuMP.@NLexpression(model, A[i, j] - B[i, j])
        end
    end
    A_minus_B
end

function Base.:*(A::Matrix{JuMP.NonlinearExpression}, B::Matrix{JuMP.NonlinearExpression})
    A_rows, A_cols = size(A)
    B_rows, B_cols = size(B)
    A_cols == B_rows || throw(DimensionMismatch("matrix sizes do not match: dimensions are $(size(A)), $(size(B))"))
    model = A[1, 1].model
    A_times_B_rows, A_times_B_cols = A_rows, B_cols
    A_times_B = Matrix{JuMP.NonlinearExpression}(undef, A_times_B_rows, A_times_B_cols)
    for j in 1:A_times_B_cols
        for i in 1:A_times_B_rows
            A_times_B[i, j] = JuMP.@NLexpression(model, sum(A[i, k] * B[k, j] for k in 1:A_cols))
        end
    end
    A_times_B
end

function checksquare(A)
    m, n = size(A)
    m == n || throw(DimensionMismatch("non square matrix of size $(size(A))"))
    return m
end

function LinearAlgebra.tr(A::Matrix{JuMP.NonlinearExpression})
    n = checksquare(A)
    model = A[1, 1].model
    JuMP.@NLexpression(model, sum(A[i, i] for i in 1:n))
end
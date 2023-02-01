# Avoid type piracy
struct NLExprMatrix
    model::JuMP.Model
    M::Matrix{JuMP.NonlinearExpression}
end

function value(A::NLExprMatrix)
    value.(A.M)
end

function *(scalar::Number, A::NLExprMatrix)
    model = A.model
    A = A.M

    m, n = size(A)
    A_scaled = Matrix{JuMP.NonlinearExpression}(undef, m, n)
    for j in 1:n
        for i in 1:m
            A_scaled[i, j] = JuMP.@NLexpression(model, scalar * A[i, j])
        end
    end

    NLExprMatrix(model, A_scaled)
end

function -(A::NLExprMatrix)
    model = A.model
    A = A.M

    m, n = size(A)
    model = A[1, 1].model
    minus_A = Matrix{JuMP.NonlinearExpression}(undef, m, n)
    for j in 1:n
        for i in 1:m
            minus_A[i, j] = JuMP.@NLexpression(model, -A[i, j])
        end
    end

    NLExprMatrix(model, minus_A)
end

function checksizematch(A, B)
    size(A) == size(B) || throw(DimensionMismatch("matrix sizes do not match: dimensions are $(size(A)), $(size(B))"))
    return size(A)
end

function +(A::NLExprMatrix, B::NLExprMatrix)
    A.model == B.model || throw(ArgumentError("matrices from two different models cannot be summed"))

    model = A.model
    A = A.M
    B = B.M

    m, n = checksizematch(A, B)
    A_plus_B = Matrix{JuMP.NonlinearExpression}(undef, m, n)
    for j in 1:n
        for i in 1:m
            A_plus_B[i, j] = JuMP.@NLexpression(model, A[i, j] + B[i, j])
        end
    end

    NLExprMatrix(model, A_plus_B)
end

function -(A::NLExprMatrix, B::NLExprMatrix)
    A.model == B.model || throw(ArgumentError("matrices from two different models cannot be summed"))

    model = A.model
    A = A.M
    B = B.M

    m, n = checksizematch(A, B)
    model = A[1, 1].model
    A_minus_B = Matrix{JuMP.NonlinearExpression}(undef, m, n)
    for j in 1:n
        for i in 1:m
            A_minus_B[i, j] = JuMP.@NLexpression(model, A[i, j] - B[i, j])
        end
    end

    NLExprMatrix(model, A_minus_B)
end

function *(A::NLExprMatrix, B::NLExprMatrix)
    A.model == B.model || throw(ArgumentError("matrices from two different models cannot be multiplied"))

    model = A.model
    A = A.M
    B = B.M

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

    NLExprMatrix(model, A_times_B)
end

function checksquare(A)
    m, n = size(A)
    m == n || throw(DimensionMismatch("non square matrix of size $(size(A))"))
    return m
end

function tr(A::NLExprMatrix)
    model = A.model
    A = A.M
    n = checksquare(A)
    JuMP.@NLexpression(model, sum(A[i, i] for i in 1:n))
end

function det(A::NLExprMatrix)
    model = A.model
    A = A.M

    n = checksquare(A)

    if n == 1
        return A[1, 1]
    elseif n == 2
        return JuMP.@NLexpression(model, A[1, 1] * A[2, 2] - A[2, 1] * A[1, 2])
    end

    minor_M = Matrix{JuMP.NonlinearExpression}(undef, (n - 1, n - 1))
    det_A = JuMP.@NLexpression(model, 0.0)

    for A_col in 1:n
        for minor_col in 1:A_col-1
            for minor_row in 1:n-1
                minor_M[minor_row, minor_col] = A[minor_row+1, minor_col]
            end
        end
        for minor_col in A_col:n-1
            for minor_row in 1:n-1
                minor_M[minor_row, minor_col] = A[minor_row+1, minor_col+1]
            end
        end

        minor = NLExprMatrix(model, minor_M)
        det_minor = det(minor)

        if isodd(A_col)
            det_A = JuMP.@NLexpression(model, det_A + A[1, A_col] * det_minor)
        else
            det_A = JuMP.@NLexpression(model, det_A - A[1, A_col] * det_minor)
        end
    end

    det_A
end

function transpose(A::NLExprMatrix)
    model = A.model
    A = A.M

    m, n = size(A)
    A_transposed = Matrix{JuMP.NonlinearExpression}(undef, (n, m))
    for j in 1:m
        for i in 1:n
            A_transposed[i, j] = A[j, i]
        end
    end

    NLExprMatrix(model, A_transposed)
end

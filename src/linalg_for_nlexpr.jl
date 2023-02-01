# Avoid type piracy
struct NLExprMatrix
    model::JuMP.Model
    internal_matrix::Matrix{JuMP.NonlinearExpression}
end

function *(scalar::Number, A::Matrix{JuMP.NonlinearExpression})
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

function -(A::Matrix{JuMP.NonlinearExpression})
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

function +(A::Matrix{JuMP.NonlinearExpression}, B::Matrix{JuMP.NonlinearExpression})
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

function -(A::Matrix{JuMP.NonlinearExpression}, B::Matrix{JuMP.NonlinearExpression})
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

function *(A::Matrix{JuMP.NonlinearExpression}, B::Matrix{JuMP.NonlinearExpression})
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

function tr(A::NLExprMatrix)
    model = A.model
    A = A.internal_matrix
    n = checksquare(A)
    JuMP.@NLexpression(model, sum(A[i, i] for i in 1:n))
end

function det(A::NLExprMatrix)
    model = A.model
    A = A.internal_matrix

    n = checksquare(A)

    if n == 1
        return A[1, 1]
    elseif n == 2
        return JuMP.@NLexpression(model, A[1, 1] * A[2, 2] - A[2, 1] * A[1, 2])
    end

    minor_internal_matrix = Matrix{JuMP.NonlinearExpression}(undef, (n - 1, n - 1))
    det_A = JuMP.@NLexpression(model, 0.0)

    for A_col in 1:n
        for minor_col in 1:A_col-1
            for minor_row in 1:n-1
                minor_internal_matrix[minor_row, minor_col] = A[minor_row+1, minor_col]
            end
        end
        for minor_col in A_col:n-1
            for minor_row in 1:n-1
                minor_internal_matrix[minor_row, minor_col] = A[minor_row+1, minor_col+1]
            end
        end

        minor = NLExprMatrix(model, minor_internal_matrix)
        det_minor = det(minor)

        if isodd(A_col)
            det_A = JuMP.@NLexpression(model, det_A + A[1, A_col] * det_minor)
        else
            det_A = JuMP.@NLexpression(model, det_A - A[1, A_col] * det_minor)
        end
    end

    det_A
end

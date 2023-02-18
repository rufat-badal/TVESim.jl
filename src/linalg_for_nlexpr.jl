# TODO: add display methods
# Avoid type piracy
struct NLExprMatrix
    model::JuMP.Model
    _matrix::Matrix{JuMP.NonlinearExpression}
end

function value(A::NLExprMatrix)
    value.(A._matrix)
end

function *(scalar::Number, A::NLExprMatrix)
    model = A.model
    A = A._matrix

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
    A = A._matrix

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
    size(A) == size(B) || throw(DimensionMismatch("sizes do not match: dimensions are $(size(A)), $(size(B))"))
    return size(A)
end

function +(A::NLExprMatrix, B::NLExprMatrix)
    A.model == B.model || throw(ArgumentError("matrices from two different models cannot be summed"))

    model = A.model
    A = A._matrix
    B = B._matrix

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
    A.model == B.model || throw(ArgumentError("matrices from two different models cannot be subtracted"))

    model = A.model
    A = A._matrix
    B = B._matrix

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
    A = A._matrix
    B = B._matrix

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
    A = A._matrix
    n = checksquare(A)
    JuMP.@NLexpression(model, sum(A[i, i] for i in 1:n))
end

function det(A::NLExprMatrix)
    n = checksquare(A._matrix)

    if n == 1
        return A._matrix[1, 1]
    elseif n == 2
        # TODO: implement getindex
        return JuMP.@NLexpression(A.model, A._matrix[1, 1] * A._matrix[2, 2] - A._matrix[2, 1] * A._matrix[1, 2])
    end

    det_A = JuMP.@NLexpression(A.model, 0.0)

    for j in 1:n
        det_minor = det(minor(A, 1, j))

        if isodd(j)
            det_A = JuMP.@NLexpression(A.model, det_A + A._matrix[1, j] * det_minor)
        else
            det_A = JuMP.@NLexpression(A.model, det_A - A._matrix[1, j] * det_minor)
        end
    end

    det_A
end

function minor(A::Matrix{JuMP.NonlinearExpression}, i, j)
    m, n = size(A)

    (1 <= i <= m && 1 <= j <= n) || throw(ArgumentError("minor indices not in the correct range, current values: ($i, $j)"))

    minor_M = Matrix{JuMP.NonlinearExpression}(undef, (m - 1, n - 1))
    if m == 1 || n == 1
        return minor_M
    end

    for minor_col in 1:j-1
        for minor_row in 1:i-1
            minor_M[minor_row, minor_col] = A[minor_row, minor_col]
        end
        for minor_row in i:m-1
            minor_M[minor_row, minor_col] = A[minor_row+1, minor_col]
        end
    end
    for minor_col in j:n-1
        for minor_row in 1:i-1
            minor_M[minor_row, minor_col] = A[minor_row, minor_col+1]
        end
        for minor_row in i:m-1
            minor_M[minor_row, minor_col] = A[minor_row+1, minor_col+1]
        end
    end

    minor_M
end

function minor(A::NLExprMatrix, i, j)
    NLExprMatrix(A.model, minor(A._matrix, i, j))
end

function transpose(A::NLExprMatrix)
    model = A.model
    A = A._matrix

    m, n = size(A)
    A_transposed = Matrix{JuMP.NonlinearExpression}(undef, (n, m))
    for j in 1:m
        for i in 1:n
            A_transposed[i, j] = A[j, i]
        end
    end

    NLExprMatrix(model, A_transposed)
end

struct NLExprVector
    model::JuMP.Model
    _vector::Vector{JuMP.NonlinearExpression}
end

function NLExprVector(model, v::Vector)
    internal_vector = Vector{JuMP.NonlinearExpression}(undef, length(v))
    for (i, x) in enumerate(v)
        internal_vector[i] = JuMP.@NLexpression(model, x)
    end

    NLExprVector(model, internal_vector)
end

function -(v::NLExprVector, w::NLExprVector)
    v.model == w.model || throw(ArgumentError("vectors from two different models cannot be summed"))
    model = v.model
    v = v._vector
    w = w._vector
    checksizematch(v, w)

    v_minus_w = [JuMP.@NLexpression(model, x - y) for (x, y) in zip(v, w)]
    NLExprVector(model, v_minus_w)
end

function +(v::NLExprVector, w::NLExprVector)
    v.model == w.model || throw(ArgumentError("vectors from two different models cannot be subtracted"))
    model = v.model
    v = v._vector
    w = w._vector
    checksizematch(v, w)

    v_plus_w = [JuMP.@NLexpression(model, x + y) for (x, y) in zip(v, w)]
    NLExprVector(model, v_plus_w)
end

function -(v::NLExprVector)
    model = v.model
    v = v._vector

    minus_v = [JuMP.@NLexpression(model, -x) for x in v]
    NLExprVector(model, minus_v)
end

function Base.hcat(v::NLExprVector, w::NLExprVector)
    v.model == w.model || throw(ArgumentError("vectors of different length cannot be concatenated: length are $(length(v)), $(length(w))"))
    model = v.model
    v = v._vector
    w = w._vector

    NLExprMatrix(model, [v w])
end
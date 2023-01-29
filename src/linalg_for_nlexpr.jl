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

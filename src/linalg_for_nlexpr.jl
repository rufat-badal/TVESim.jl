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

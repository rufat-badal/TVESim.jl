struct JuMPExpression
    model::JuMP.Model
    expr
end

Base.show(io::IO, x::JuMPExpression) = show(io, x.expr)

jumpexpression_array(model, X::Array) = [JuMPExpression(model, x) for x in X]

JuMP.value(x::JuMPExpression) = JuMP.value(JuMP.add_nonlinear_expression(x.model, x.expr))

JuMPExpression(x::Number) = x

Base.zero(::Type{JuMPExpression}) = 0 # needed for LinearAlgebra.tr
Base.zero(x::JuMPExpression) = JuMPExpression(x.model, 0) # needed for LinearAlgebra.dot

Base.length(x::JuMPExpression) = 1

Base.iterate(x::JuMPExpression) = (x, nothing)
Base.iterate(x::JuMPExpression, ::Any) = nothing

Base.transpose(x::JuMPExpression) = x

Base.:^(x::JuMPExpression, power::Integer) = JuMPExpression(x.model, :($(x.expr)^$(power)))

Base.:-(x::JuMPExpression) = JuMPExpression(x.model, :(-$(x.expr)))
Base.inv(x::JuMPExpression) = sJuMPExpression(x.model, :(1 / $(x.expr)))

Base.:+(x::JuMPExpression, y::JuMPExpression) = JuMPExpression(x.model, :($(x.expr) + $(y.expr)))
Base.:+(x::JuMPExpression, y) = JuMPExpression(x.model, :($(x.expr) + $(y)))
Base.:+(x, y::JuMPExpression) = JuMPExpression(y.model, :($(x) + $(y.expr)))

Base.:-(x::JuMPExpression, y::JuMPExpression) = JuMPExpression(x.model, :($(x.expr) - $(y.expr)))
Base.:-(x::JuMPExpression, y) = JuMPExpression(x.model, :($(x.expr) - $(y)))
Base.:-(x, y::JuMPExpression) = JuMPExpression(y.model, :($(x) - $(y.expr)))

Base.:*(x::JuMPExpression, y::JuMPExpression) = JuMPExpression(x.model, :($(x.expr) * $(y.expr)))
Base.:*(x::JuMPExpression, y) = JuMPExpression(x.model, :($(x.expr) * $(y)))
Base.:*(x, y::JuMPExpression) = JuMPExpression(y.model, :($(x) * $(y.expr)))

Base.:/(x::JuMPExpression, y::JuMPExpression) = JuMPExpression(x.model, :($(x.expr) / $(y.expr)))
Base.:/(x::JuMPExpression, y) = JuMPExpression(x.model, :($(x.expr) / $(y)))
Base.:/(x, y::JuMPExpression) = JuMPExpression(y.model, :($(x) / $(y.expr)))

Base.:*(λ::JuMPExpression, X::Matrix{JuMPExpression}) = [λ * x for x in X]
Base.:*(λ::JuMP.AbstractJuMPScalar, X::Matrix{JuMPExpression}) = [λ * x for x in X]
Base.:*(X::Matrix{JuMPExpression}, λ::JuMPExpression) = [x * λ for x in X]
Base.:*(X::Matrix{JuMPExpression}, λ::JuMP.AbstractJuMPScalar) = [x * λ for x in X]

Base.:/(X::Matrix{JuMPExpression}, λ::JuMPExpression) = [x / λ for x in X]
Base.:/(X::Matrix{JuMPExpression}, λ::JuMP.AbstractJuMPScalar) = [x / λ for x in X]

LinearAlgebra.dot(x::JuMPExpression, y::JuMPExpression) = x * y
LinearAlgebra.dot(x, y::JuMPExpression) = x * y
LinearAlgebra.dot(x::JuMPExpression, y) = x * y

function checksquare(A)
    m, n = size(A)
    m == n || throw(DimensionMismatch("operation does not support non-square matrices, matrix has size $(size(A))"))
    return m
end

function minor(X::Matrix{T}, i, j) where T
    m, n = size(X)

    (1 <= i <= m && 1 <= j <= n) || throw(ArgumentError("minor indices not in the correct range, current values: ($i, $j)"))

    Xminor = Matrix{T}(undef, m - 1, n - 1)
    if m == 1 || n == 1
        return Xminor
    end

    for Xminor_col in 1:j-1
        for Xminor_row in 1:i-1
            Xminor[Xminor_row, Xminor_col] = X[Xminor_row, Xminor_col]
        end
        for Xminor_row in i:m-1
            Xminor[Xminor_row, Xminor_col] = X[Xminor_row+1, Xminor_col]
        end
    end
    for Xminor_col in j:n-1
        for Xminor_row in 1:i-1
            Xminor[Xminor_row, Xminor_col] = X[Xminor_row, Xminor_col+1]
        end
        for Xminor_row in i:m-1
            Xminor[Xminor_row, Xminor_col] = X[Xminor_row+1, Xminor_col+1]
        end
    end

    Xminor
end

function det(X::Matrix{JuMPExpression})
    n = checksquare(X)

    if n == 1
        return X[1, 1]
    elseif n == 2
        return X[1, 1] * X[2, 2] - X[2, 1] * X[1, 2]
    end

    det_X = zero(JuMPExpression)

    for j in 1:n
        det_minor = det(minor(X, 1, j))

        if isodd(j)
            det_X = det_X + X[1, j] * det_minor
        else
            det_X = det_X - X[1, j] * det_minor
        end
    end

    det_X
end
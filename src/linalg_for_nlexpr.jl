struct AdvancedNonlinearExpression
    model::JuMP.Model
    _expression::JuMP.NonlinearExpression
end

function AdvancedNonlinearExpression(model::JuMP.Model, val)
    AdvancedNonlinearExpression(
        model,
        JuMP.@NLexpression(model, val)
    )
end

function Base.show(io::IO, x::AdvancedNonlinearExpression)
    show(io, x._expression)
end

function Base.:-(x::AdvancedNonlinearExpression)
    AdvancedNonlinearExpression(x.model, JuMP.@NLexpression(x.model, -x._expression))
end

function checkmodelmatch(x, y)
    x.model == y.model || throw(ArgumentError("nonlinear expressions belong to different models"))
end

function Base.:+(x::AdvancedNonlinearExpression, y::AdvancedNonlinearExpression)
    checkmodelmatch(x, y)
    AdvancedNonlinearExpression(x.model, JuMP.@NLexpression(x.model, x._expression + y._expression))
end

function Base.:+(x::AdvancedNonlinearExpression, y::Number)
    AdvancedNonlinearExpression(x.model, JuMP.@NLexpression(x.model, x._expression + y))
end

function Base.:+(x::Number, y::AdvancedNonlinearExpression)
    y + x
end

function Base.:-(x::AdvancedNonlinearExpression, y::AdvancedNonlinearExpression)
    checkmodelmatch(x, y)
    AdvancedNonlinearExpression(x.model, JuMP.@NLexpression(x.model, x._expression - y._expression))
end

function Base.:-(x::AdvancedNonlinearExpression, y::Number)
    AdvancedNonlinearExpression(x.model, JuMP.@NLexpression(x.model, x._expression - y))
end

function Base.:-(x::Number, y::AdvancedNonlinearExpression)
    AdvancedNonlinearExpression(y.model, JuMP.@NLexpression(y.model, x - y._expression))
end

function Base.:*(x::AdvancedNonlinearExpression, y::AdvancedNonlinearExpression)
    checkmodelmatch(x, y)
    AdvancedNonlinearExpression(x.model, JuMP.@NLexpression(x.model, x._expression * y._expression))
end

function Base.:*(x::AdvancedNonlinearExpression, y::Number)
    AdvancedNonlinearExpression(x.model, JuMP.@NLexpression(x.model, y * x._expression))
end

function Base.:*(x::Number, y::AdvancedNonlinearExpression)
    y * x
end

transpose(x::AdvancedNonlinearExpression) = x

function Base.:^(x::AdvancedNonlinearExpression, power::Int)
    AdvancedNonlinearExpression(x.model, JuMP.@NLexpression(x.model, x._expression^power))
end

function dot(X::Matrix{AdvancedNonlinearExpression}, Y::Matrix{AdvancedNonlinearExpression})
    model = X[1, 1].model
    AdvancedNonlinearExpression(
        model,
        JuMP.@NLexpression(
            model,
            sum(x._expression * y._expression for (x, y) in zip(X, Y))
        )
    )
end

function dot(X::Matrix{AdvancedNonlinearExpression}, Y::Matrix)
    model = X[1, 1].model
    AdvancedNonlinearExpression(
        model,
        JuMP.@NLexpression(
            model,
            sum(x._expression * y for (x, y) in zip(X, Y))
        )
    )
end

function dot(X::Matrix, Y::Matrix{AdvancedNonlinearExpression})
    dot(Y, X)
end

function norm_sqr(X::Matrix{AdvancedNonlinearExpression})
    model = X[1, 1].model
    AdvancedNonlinearExpression(
        model,
        JuMP.@NLexpression(
            model,
            sum(x._expression^2 for x in X)
        )
    )
end

function get_product_size(A, B)
    
end

function Base.:*(X::Matrix{AdvancedNonlinearExpression}, Y::Matrix{AdvancedNonlinearExpression})
    X_rows, X_cols = size(X)
    Y_rows, Y_cols = size(Y)
    X_cols == Y_rows || throw(DimensionMismatch("matrix sizes do not match: dimensions are $(size(A)), $(size(B))"))
    XY_rows, XY_cols = X_rows, Y_cols

    model = X[1, 1].model

    XY = Matrix{AdvancedNonlinearExpression}(undef, XY_rows, XY_cols)

    for j in 1:XY_cols
        for i in 1:XY_rows
            XY[i, j] = AdvancedNonlinearExpression(
                model,
                JuMP.@NLexpression(model, sum(X[i, k]._expression * Y[k, j]._expression for k in 1:X_cols))
            )
        end
    end

    XY
end

function Base.:*(X::Matrix{AdvancedNonlinearExpression}, Y::Matrix)
    X_rows, X_cols = size(X)
    Y_rows, Y_cols = size(Y)
    X_cols == Y_rows || throw(DimensionMismatch("matrix sizes do not match: dimensions are $(size(A)), $(size(B))"))
    XY_rows, XY_cols = X_rows, Y_cols

    model = X[1, 1].model

    XY = Matrix{AdvancedNonlinearExpression}(undef, XY_rows, XY_cols)

    for j in 1:XY_cols
        for i in 1:XY_rows
            XY[i, j] = AdvancedNonlinearExpression(
                model,
                JuMP.@NLexpression(model, sum(X[i, k]._expression * Y[k, j] for k in 1:X_cols))
            )
        end
    end

    XY
end

function Base.:*(X::Matrix, Y::Matrix{AdvancedNonlinearExpression})
    X_rows, X_cols = size(X)
    Y_rows, Y_cols = size(Y)
    X_cols == Y_rows || throw(DimensionMismatch("matrix sizes do not match: dimensions are $(size(A)), $(size(B))"))
    XY_rows, XY_cols = X_rows, Y_cols

    model = Y[1, 1].model

    XY = Matrix{AdvancedNonlinearExpression}(undef, XY_rows, XY_cols)

    for j in 1:XY_cols
        for i in 1:XY_rows
            XY[i, j] = AdvancedNonlinearExpression(
                model,
                JuMP.@NLexpression(model, sum(X[i, k] * Y[k, j]._expression for k in 1:X_cols))
            )
        end
    end

    XY
end

function checksquare(A)
    m, n = size(A)
    m == n || throw(DimensionMismatch("operation does not support non-square matrices, matrix has size $(size(A))"))
    return m
end

function tr(X::Matrix{AdvancedNonlinearExpression})
    model = X[1, 1].model
    n = checksquare(X)
    AdvancedNonlinearExpression(
        model,
        JuMP.@NLexpression(model, sum(X[i, i]._expression for i in 1:n))
    )
end

function minor(X::Matrix{AdvancedNonlinearExpression}, i, j)
    m, n = size(X)

    (1 <= i <= m && 1 <= j <= n) || throw(ArgumentError("minor indices not in the correct range, current values: ($i, $j)"))

    Xminor = Matrix{AdvancedNonlinearExpression}(undef, m - 1, n - 1)
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

function det(X::Matrix{AdvancedNonlinearExpression})
    n = checksquare(X)
    model = X[1, 1].model

    if n == 1
        return X[1, 1]
    elseif n == 2
        return AdvancedNonlinearExpression(
            model,
            JuMP.@NLexpression(model, X[1, 1]._expression * X[2, 2]._expression - X[2, 1]._expression * X[1, 2]._expression)
        )
    end

    det_X = AdvancedNonlinearExpression(model, 0)

    for j in 1:n
        det_minor = det(minor(X, 1, j))

        if isodd(j)
            det_X = AdvancedNonlinearExpression(
                model,
                JuMP.@NLexpression(model, det_X._expression + X[1, j]._expression * det_minor._expression)
            )
        else
            det_X = AdvancedNonlinearExpression(
                model,
                JuMP.@NLexpression(model, det_X._expression - X[1, j]._expression * det_minor._expression)
            )
        end
    end

    det_X
end

function adjugate(X::Matrix{AdvancedNonlinearExpression})
    n = checksquare(X)
    Xadj = Matrix{AdvancedNonlinearExpression}(undef, n, n)

    for j in 1:n
        for i in 1:n
            det_Xminor = det(minor(X, j, i))
            if iseven(i + j)
                Xadj[i, j] = det_Xminor
            else
                Xadj[i, j] = -det_Xminor
            end
        end
    end

    Xadj
end


struct AdvancedNonlinearExpression
    model::JuMP.Model
    _expression::JuMP.NonlinearExpression
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

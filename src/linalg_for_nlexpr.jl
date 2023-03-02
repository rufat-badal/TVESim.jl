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
    y - x
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

function dot(X::AbstractMatrix{AdvancedNonlinearExpression}, Y::AbstractMatrix{AdvancedNonlinearExpression})
    model = X[1, 1].model
    AdvancedNonlinearExpression(
        model,
        JuMP.@NLexpression(
            model,
            sum(x._expression * y._expression for (x, y) in zip(X, Y))
        )
    )
end

function norm_sqr(X::AbstractMatrix{AdvancedNonlinearExpression})
    model = X[1, 1].model
    AdvancedNonlinearExpression(
        model,
        JuMP.@NLexpression(
            model,
            sum(x._expression^2 for x in X)
        )
    )
end
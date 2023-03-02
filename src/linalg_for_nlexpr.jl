struct AdvancedNonlinearExpression
    model::JuMP.Model
    _expression::JuMP.NonlinearExpression
end

function Base.show(io::IO, e::AdvancedNonlinearExpression)
    show(io, e._expression)
end

function Base.:-(e::AdvancedNonlinearExpression)
    AdvancedNonlinearExpression(e.model, JuMP.@NLexpression(e.model, -e._expression))
end

function checkmodelmatch(e1, e2)
    e1.model == e2.model || throw(ArgumentError("nonlinear expressions belong to different models"))
end

function Base.:+(e1::AdvancedNonlinearExpression, e2::AdvancedNonlinearExpression)
    checkmodelmatch(e1, e2)
    AdvancedNonlinearExpression(e1.model, JuMP.@NLexpression(e1.model, e1._expression + e2._expression))
end

function Base.:+(e1::AdvancedNonlinearExpression, e2::Number)
    AdvancedNonlinearExpression(e1.model, JuMP.@NLexpression(e1.model, e1._expression + e2))
end

function Base.:+(e1::Number, e2::AdvancedNonlinearExpression)
    e2 + e1
end

function Base.:-(e1::AdvancedNonlinearExpression, e2::AdvancedNonlinearExpression)
    checkmodelmatch(e1, e2)
    AdvancedNonlinearExpression(e1.model, JuMP.@NLexpression(e1.model, e1._expression - e2._expression))
end

function Base.:-(e1::AdvancedNonlinearExpression, e2::Number)
    AdvancedNonlinearExpression(e1.model, JuMP.@NLexpression(e1.model, e1._expression - e2))
end

function Base.:-(e1::Number, e2::AdvancedNonlinearExpression)
    e2 - e1
end

function Base.:*(e1::AdvancedNonlinearExpression, e2::AdvancedNonlinearExpression)
    checkmodelmatch(e1, e2)
    AdvancedNonlinearExpression(e1.model, JuMP.@NLexpression(e1.model, e1._expression * e2._expression))
end

function Base.:*(e1::AdvancedNonlinearExpression, e2::Number)
    AdvancedNonlinearExpression(e1.model, JuMP.@NLexpression(e1.model, e2 * e1._expression))
end

function Base.:*(e1::Number, e2::AdvancedNonlinearExpression)
    e2 * e1
end

transpose(e::AdvancedNonlinearExpression) = e

function Base.:^(e::AdvancedNonlinearExpression, power::Int)
    AdvancedNonlinearExpression(e.model, JuMP.@NLexpression(e.model, e._expression ^ power))
end
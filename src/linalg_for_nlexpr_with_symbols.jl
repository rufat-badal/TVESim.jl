struct JuMPExpression
    model::JuMP.Model
    expr 
end

function Base.show(io::IO, x::JuMPExpression)
    show(io, x.expr)
end

function jumpexpression_array(model, X::Array)
    Xout = [JuMPExpression(model, x) for x in X]
end

function JuMP.add_nonlinear_expression(model::JuMP.Model, x::JuMPExpression)
    JuMP.add_nonlinear_expression(model, x.expr)
end

Base.length(x::JuMPExpression) = 1

Base.iterate(x::JuMPExpression) = (x, nothing)
Base.iterate(x::JuMPExpression, ::Any) = nothing

function Base.:-(x::JuMPExpression)
    JuMPExpression(x.model, :(-$(x.expr)))
end

function Base.:+(x::JuMPExpression, y::JuMPExpression)
    JuMPExpression(x.model, :($(x.expr) + $(y.expr)))
end

function Base.:+(x::JuMPExpression, y)
    JuMPExpression(x.model, :($(x.expr) + $(y)))
end

function Base.:+(x, y::JuMPExpression)
    JuMPExpression(y.model, :($(x) + $(y.expr)))
end

function Base.:-(x::JuMPExpression, y::JuMPExpression)
    JuMPExpression(x.model, :($(x.expr) - $(y.expr)))
end

function Base.:-(x::JuMPExpression, y)
    JuMPExpression(x.model, :($(x.expr) - $(y)))
end

function Base.:-(x, y::JuMPExpression)
    JuMPExpression(y.model, :($(x) - $(y.expr)))
end

function Base.:*(x::JuMPExpression, y::JuMPExpression)
    JuMPExpression(x.model, :($(x.expr) * $(y.expr)))
end

function Base.:*(x::JuMPExpression, y)
    JuMPExpression(x.model, :($(x.expr) * $(y)))
end

function Base.:*(x, y::JuMPExpression)
    JuMPExpression(y.model, :($(x) * $(y.expr)))
end

function Base.:*(λ::JuMPExpression, X::Matrix{JuMPExpression})
    λ .* X
end

function Base.:*(λ::JuMP.AbstractJuMPScalar, X::Matrix{JuMPExpression})
    [λ * x for x in X]
end

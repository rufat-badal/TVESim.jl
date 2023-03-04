struct JuMPExpression
    model::JuMP.Model
    expr 
end

function Base.show(io::IO, x::JuMPExpression)
    show(io, x.expr)
end

function jumpexpression_array(model, X::Array)
    [JuMPExpression(model, x) for x in X]
end

function JuMP.value(x::JuMPExpression)
    JuMP.value(JuMP.add_nonlinear_expression(x.model, x.expr))
end

JuMPExpression(x::Number) = x

Base.zero(::Type{JuMPExpression}) = 0 # needed for LinearAlgebra.tr
Base.zero(x::JuMPExpression) = JuMPExpression(x.model, 0) # needed for LinearAlgebra.dot

Base.length(x::JuMPExpression) = 1

Base.iterate(x::JuMPExpression) = (x, nothing)
Base.iterate(x::JuMPExpression, ::Any) = nothing

Base.transpose(x::JuMPExpression) = x

function Base.:^(x::JuMPExpression, power::Integer)
    JuMPExpression(x.model, :($(x.expr)^$(power)))
end

function Base.:-(x::JuMPExpression)
    JuMPExpression(x.model, :(-$(x.expr)))
end

function Base.inv(x::JuMPExpression)
    JuMPExpression(x.model, :(1 / $(x.expr)))
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

function Base.:/(x::JuMPExpression, y::JuMPExpression)
    JuMPExpression(x.model, :($(x.expr) / $(y.expr)))
end

function Base.:/(x::JuMPExpression, y)
    JuMPExpression(x.model, :($(x.expr) / $(y)))
end

function Base.:/(x, y::JuMPExpression)
    JuMPExpression(y.model, :($(x) / $(y.expr)))
end

function Base.:*(λ::JuMPExpression, X::Matrix{JuMPExpression})
    [λ * x for x in X]
end

function Base.:*(λ::JuMP.AbstractJuMPScalar, X::Matrix{JuMPExpression})
    [λ * x for x in X]
end

function Base.:*(X::Matrix{JuMPExpression}, λ::JuMPExpression)
   [x * λ for x in X] 
end

function Base.:*(X::Matrix{JuMPExpression}, λ::JuMP.AbstractJuMPScalar)
   [x * λ for x in X] 
end

function Base.:/(X::Matrix{JuMPExpression}, λ::JuMPExpression)
   [x / λ for x in X] 
end

function Base.:/(X::Matrix{JuMPExpression}, λ::JuMP.AbstractJuMPScalar)
   [x / λ for x in X] 
end

function LinearAlgebra.dot(x::JuMPExpression, y::JuMPExpression)
    x * y
end

function LinearAlgebra.dot(x, y::JuMPExpression)
    x * y
end

function LinearAlgebra.dot(x::JuMPExpression, y)
    x * y
end
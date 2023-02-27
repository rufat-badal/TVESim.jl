function austenite_percentage(θ)
    return 1 - 1 / (1 + θ)
end

function neo_hook(F::NLExprMatrix)
    trace_C = tr(transpose(F) * F)
    det_F = det(F)
    JuMP.@NLexpression(F.model, trace_C - 2 - 2 * log(det_F) + (det_F - 1)^2)
end

function gradient_austenite_potential(F::NLExprMatrix)
    det_F = det(F)
    scalar_expr = JuMP.@NLexpression(F.model, det_F - 1 / det_F - 1)
    det_gradient = NLExprMatrix([F[2, 2] -F[2, 1]; -F[1, 2] F[1, 1]])
    2 * (F + scalar_expr * det_gradient)
end

function gradient_martensite_potential(F::NLExprMatrix, scaling_matrix::Matrix)
    # chain rule
    gradient_austenite_potential(F * scaling_matrix) * transpose(scaling_matrix)
end
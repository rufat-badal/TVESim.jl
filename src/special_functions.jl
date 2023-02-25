function austenite_percentage(θ)
    return 1 - 1 / (1 + θ)
end

function neo_hook(F::NLExprMatrix)
    trace_C = tr(transpose(F) * F)
    det_F = det(F)
    JuMP.@NLexpression(F.model, trace_C - 2 - 2 * log(det_F) + (det_F - 1)^2)
end
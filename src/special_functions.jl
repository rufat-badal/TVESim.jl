austenite_percentage(θ) = 1 - 1 / (1 + θ)

dissipation_potential(dot_C::Matrix{JuMPExpression}) = norm_sqr(dot_C) / 2

function neo_hook(F::Matrix{JuMPExpression})
    trace_C = tr(transpose(F) * F)
    det_F = det(F)
    trace_C - 2 - 2 * log(det_F) + (det_F - 1)^2
end

function gradient_austenite_potential(F::Matrix{JuMPExpression})
    det_F = det(F)
    gradient_det = [F[2, 2] -F[2, 1]; -F[1, 2] F[1, 1]]
    2 * (F + (det_F - 1 / det_F - 1) * gradient_det)
end

function gradient_martensite_potential(F::Matrix{JuMPExpression}, scaling_matrix::Matrix)
    # chain rule
    gradient_austenite_potential(F * scaling_matrix) * transpose(scaling_matrix)
end

austenite_potential(F::Matrix{JuMPExpression}) = neo_hook(F)
martensite_potential(F, scaling_matrix) = austenite_potential(F * scaling_matrix)

function internal_energy_weight(θ)
    # a = austenite_percentage
    # a(θ) - θ a'(θ)
    θ^2 / (1 + θ)^2
end

function antider_internal_energy_weight(θ)
    # a = austenite_percentage
    # integral of a(s) - s a'(s) from s = 0 to s = θ
    (θ * (2 + θ)) / (1 + θ) - 2 * log(1 + θ)
end

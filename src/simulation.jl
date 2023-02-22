struct SimulationStep
    x::Vector{Float64}
    y::Vector{Float64}
    θ::Vector{Float64}
end

struct MechanicalStep
    model::JuMP.Model
    prev_x::Vector{JuMP.NonlinearParameter}
    prev_y::Vector{JuMP.NonlinearParameter}
    prev_θ::Vector{JuMP.NonlinearParameter}
    x::Vector{JuMP.VariableRef}
    y::Vector{JuMP.VariableRef}
end

function MechanicalStep(grid::SimulationGrid, search_rad)
    m = JuMP.Model(() -> MadNLP.Optimizer(print_level=MadNLP.WARN, blas_num_threads=8))

    # previous steps
    num_vertices = grid.num_vertices
    JuMP.@NLparameter(m, prev_x[i=1:num_vertices] == grid.x[i])
    JuMP.@NLparameter(m, prev_y[i=1:num_vertices] == grid.y[i])
    JuMP.@NLparameter(m, prev_θ[i=1:num_vertices] == grid.θ[i])

    # variables
    JuMP.@variable(
        m,
        value(prev_x[i]) - search_rad <= x[i=1:num_vertices] <= value(prev_x[i]) + search_rad,
        start = value(prev_x[i])
    )
    JuMP.@constraint(m, fix_x[i=1:grid.num_dirichlet_vertices], x[i] == grid.x[i])
    JuMP.@variable(
        m,
        value(prev_y[i]) - search_rad <= y[i=1:num_vertices] <= value(prev_y[i]) + search_rad,
        start = value(prev_y[i])
    )
    JuMP.@constraint(m, fix_y[i=1:grid.num_dirichlet_vertices], y[i] == grid.y[i])

    MechanicalStep(m, prev_x, prev_y, prev_θ, x, y)
end

Base.@kwdef struct Simulation
    grid::SimulationGrid
    deformation_search_radius::Number
    temperature_search_radius::Number
    fps::Number = 30
    initial_temperature::Number = 0.0
    initial_scaling::Number = 1.0
    shape_memory_scaling::Number = 2.0
    heat_conductivity::Vector{Number} = [1.0, 1.0]
    heat_transfer_coefficient::Number = 0.5
    entropic_heat_capacity::Number = 10.0
    external_temperature::Number = 0.0
    mechanical_step::MechanicalStep
    steps::Vector{SimulationStep}
end

function Simulation(grid::SimulationGrid, deformation_search_radius, temperature_search_radius)
    mechanical_step = MechanicalStep(grid, deformation_search_radius)
    steps = [SimulationStep(value.(mechanical_step.prev_x), value.(mechanical_step.prev_y), value.(mechanical_step.prev_θ))]
    simulation = Simulation(; grid, deformation_search_radius, temperature_search_radius, mechanical_step, steps)
    run_first_mechanical_step(simulation)
    push!(steps, SimulationStep(value.(mechanical_step.x), value.(mechanical_step.y), value.(mechanical_step.prev_θ)))
    plot(steps[end], grid.triangles)

    # simulation
end

function plot(step::SimulationStep, triangles)
    plot_width = 1000
    strokewidth = 1
    padding_perc = 0.5

    x, y = step.x, step.y
    vertex_coords = [x y]
    faces = Matrix{Int}(undef, length(triangles), 3)
    for (i, T) in enumerate(triangles)
        for (j, v) in enumerate(T)
            faces[i, j] = v
        end
    end

    min_x, max_x = minimum(x), maximum(x)
    min_y, max_y = minimum(y), maximum(y)
    width = max_x - min_x
    height = max_y - min_y
    aspect = width / height
    plot_height = plot_width / aspect
    fig = CairoMakie.Figure(resolution=(plot_width, plot_height))
    horizontal_padding = padding_perc / 100 * width
    vertical_padding = padding_perc / 100 * height
    ax = CairoMakie.Axis(
        fig[1, 1],
        limits=(
            min_x - horizontal_padding, max_x + vertical_padding,
            min_y - horizontal_padding, max_y + horizontal_padding),
        aspect=aspect)
    CairoMakie.hidedecorations!(ax)
    CairoMakie.hidespines!(ax)
    CairoMakie.poly!(vertex_coords, faces, color=:transparent, strokewidth=strokewidth, shading=true)

    fig
end

function run_first_mechanical_step(simulation::Simulation)
    create_mechanical_step_objective(simulation)
    JuMP.optimize!(simulation.mechanical_step.model)
end

function create_mechanical_step_objective(simulation::Simulation)
    m = simulation.mechanical_step.model
    prev_x = simulation.mechanical_step.prev_x
    prev_y = simulation.mechanical_step.prev_y
    prev_θ = simulation.mechanical_step.prev_θ
    x = simulation.mechanical_step.x
    y = simulation.mechanical_step.y
    grid = simulation.grid

    # compute strains and symmetrized strain-rates
    prev_triangles = [[NLExprVector(m, [prev_x[i], prev_y[i]]) for i in T] for T in grid.triangles]
    triangles = [[NLExprVector(m, [x[i], y[i]]) for i in T] for T in grid.triangles]
    reference_triangles = [[[grid.x[i], grid.y[i]] for i in T] for T in grid.triangles]
    prev_strains = strain.(zip(prev_triangles, reference_triangles))
    strains = strain.(zip(triangles, reference_triangles))
    strain_rates = [F - prev_F for (F, prev_F) in zip(strains, prev_strains)]
    symmetrized_strain_rates = [
        transpose(dot_F) * prev_F + transpose(prev_F) * dot_F
        for (prev_F, dot_F) in zip(prev_strains, strain_rates)
    ]
    JuMP.register(m, :austenite_percentage, 1, austenite_percentage; autodiff=true)
    austenite_percentages = Vector{JuMP.NonlinearExpression}(undef, length(grid.triangles))
    for (i, (i1, i2, i3)) in enumerate(grid.triangles)
        austenite_percentages[i] = JuMP.@NLexpression(
            m,
            (austenite_percentage(prev_θ[i1]) + austenite_percentage(prev_θ[i2]) + austenite_percentage(prev_θ[i3])) / 3
        )
    end

    # objective
    scaling_matrix = [1/simulation.shape_memory_scaling 0; 0 1]
    scaled_strains = [F * scaling_matrix for F in strains]
    JuMP.@NLexpression(
        m, elastic_energy,
        0.5 * sum(
            (1 - a_perc) * neo_hook_F_scaled + a_perc * neo_hook_F
            for (F, a_perc, neo_hook_F, neo_hook_F_scaled) in zip(
                strains,
                austenite_percentages,
                neo_hook.(strains),
                neo_hook.(scaled_strains)
            )
        )
    )
    JuMP.@NLexpression(m, dissipation, 0.5 * sum(d for d in norm_sqr.(symmetrized_strain_rates)))
    JuMP.@NLobjective(m, Min, elastic_energy + simulation.fps * dissipation)
end

function neo_hook(F::NLExprMatrix)
    trace_C = tr(transpose(F) * F)
    det_F = det(F)
    JuMP.@NLexpression(F.model, trace_C - 2 - 2 * log(det_F) + (det_F - 1)^2)
end

function austenite_percentage(θ)
    return 1 - 1 / (1 + θ)
end

function strain((triangle, reference_triangle))
    a, b, c = reference_triangle
    gradient_equilateral_to_reference = [b - a c - a]
    a, b, c = triangle
    gradient_equilateral_to_current = [b - a c - a]

    gradient_equilateral_to_current * inv(gradient_equilateral_to_reference)
end
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

function MechanicalStep(grid::SimulationGrid, search_rad::Number)
    m = JuMP.Model(() -> MadNLP.Optimizer(print_level=MadNLP.WARN, blas_num_threads=16))

    # previous steps
    num_vertices = grid.num_vertices
    JuMP.@NLparameter(m, prev_x[i=1:num_vertices] == grid.x[i])
    JuMP.@NLparameter(m, prev_y[i=1:num_vertices] == grid.y[i])
    JuMP.@NLparameter(m, prev_θ[i=1:num_vertices] == grid.θ[i])

    # variables
    JuMP.@variable(
        m,
        JuMP.value(prev_x[i]) - search_rad <= x[i=1:num_vertices] <= JuMP.value(prev_x[i]) + search_rad,
        start = JuMP.value(prev_x[i])
    )
    JuMP.@constraint(m, fix_x[i=1:grid.num_dirichlet_vertices], x[i] == grid.x[i])
    JuMP.@variable(
        m,
        JuMP.value(prev_y[i]) - search_rad <= y[i=1:num_vertices] <= JuMP.value(prev_y[i]) + search_rad,
        start = JuMP.value(prev_y[i])
    )
    JuMP.@constraint(m, fix_y[i=1:grid.num_dirichlet_vertices], y[i] == grid.y[i])

    MechanicalStep(m, prev_x, prev_y, prev_θ, x, y)
end

struct ThermalStep
    model::JuMP.Model
    prev_x::Vector{JuMP.NonlinearParameter}
    prev_y::Vector{JuMP.NonlinearParameter}
    prev_θ::Vector{JuMP.NonlinearParameter}
    x::Vector{JuMP.NonlinearParameter}
    y::Vector{JuMP.NonlinearParameter}
    θ::Vector{JuMP.VariableRef}
end

function ThermalStep(grid::SimulationGrid, mechanical_step::MechanicalStep, search_rad::Number)
    m = JuMP.Model(() -> MadNLP.Optimizer(print_level=MadNLP.WARN, blas_num_threads=8))

    num_vertices = grid.num_vertices
    JuMP.@NLparameter(m, prev_x[i=1:num_vertices] == JuMP.value(mechanical_step.prev_x[i]))
    JuMP.@NLparameter(m, prev_y[i=1:num_vertices] == JuMP.value(mechanical_step.prev_y[i]))
    JuMP.@NLparameter(m, prev_θ[i=1:num_vertices] == JuMP.value(mechanical_step.prev_θ[i]))
    JuMP.@NLparameter(m, x[i=1:num_vertices] == JuMP.value(mechanical_step.x[i]))
    JuMP.@NLparameter(m, y[i=1:num_vertices] == JuMP.value(mechanical_step.y[i]))

    JuMP.@variable(
        m,
        JuMP.value(prev_θ[i]) - search_rad <= θ[i=1:num_vertices] <= JuMP.value(prev_θ[i]) + search_rad,
        start = JuMP.value(prev_θ[i])
    )

    ThermalStep(m, prev_x, prev_y, prev_θ, x, y, θ)
end

struct Simulation
    grid::SimulationGrid
    deformation_search_radius::Number
    temperature_search_radius::Number
    shape_memory_scaling::Number
    initial_temperature::Number
    fps::Number
    heat_transfer_coefficient::Number
    heat_conductivity::Vector{Number}
    entropic_heat_capacity::Number
    external_temperature::Number
    mechanical_step::MechanicalStep
    steps::Vector{SimulationStep}
end

function Simulation(
    grid::SimulationGrid;
    shape_memory_scaling=1.5,
    initial_temperature=0,
    fps=30,
    heat_transfer_coefficient=0.5,
    heat_conductivity=[1.0, 1.0],
    entropic_heat_capacity=10.0,
    external_temperature=0.0
)
    width = maximum(grid.x) - minimum(grid.x)
    height = maximum(grid.y) - minimum(grid.y)
    deformation_search_radius = 1.1 * shape_memory_scaling * max(width, height)

    mechanical_step = MechanicalStep(grid, deformation_search_radius)
    x = JuMP.value.(mechanical_step.prev_x)
    y = JuMP.value.(mechanical_step.prev_y)
    θ = JuMP.value.(mechanical_step.prev_θ)
    steps = [SimulationStep(x, y, θ)]
    create_objective!(mechanical_step, grid, shape_memory_scaling, fps)
    JuMP.optimize!(mechanical_step.model)

    temperature_search_radius = initial_temperature + 10
    thermal_step = ThermalStep(grid, mechanical_step, temperature_search_radius)
    create_objective!(thermal_step, grid, heat_transfer_coefficient, heat_conductivity, entropic_heat_capacity, external_temperature)

    x = JuMP.value.(mechanical_step.x)
    y = JuMP.value.(mechanical_step.y)
    θ = JuMP.value.(mechanical_step.prev_θ)
    push!(steps, SimulationStep(x, y, θ))

    Simulation(
        grid,
        deformation_search_radius,
        temperature_search_radius,
        shape_memory_scaling,
        initial_temperature,
        fps,
        heat_transfer_coefficient,
        heat_conductivity,
        entropic_heat_capacity,
        external_temperature,
        mechanical_step,
        steps
    )
end

function create_objective!(
    thermal_step::ThermalStep, grid::SimulationGrid,
    heat_transfer_coefficient, heat_conductivity, entropic_heat_capacity, external_temperature
)
    m = thermal_step.model
    prev_x = thermal_step.prev_x
    prev_y = thermal_step.prev_y
    prev_θ = thermal_step.prev_θ
    x = thermal_step.x
    y = thermal_step.y
    θ = thermal_step.θ

    # recompute strain and strain rate
    prev_strains, strains = get_strains(m, prev_x, prev_y, x, y, grid)
    symmetrized_strain_rates = get_symmetrized_strain_rates(prev_strains, strains)
end

function update!(mechanical_step::MechanicalStep)
    JuMP.set_value.(mechanical_step.prev_x, JuMP.value.(mechanical_step.x))
    JuMP.set_value.(mechanical_step.prev_y, JuMP.value.(mechanical_step.y))
end

function solve!(mechanical_step::MechanicalStep)
    JuMP.optimize!(mechanical_step.model)
end

function append_step!(simulation::Simulation)
    x = JuMP.value.(simulation.mechanical_step.x)
    y = JuMP.value.(simulation.mechanical_step.y)
    θ = JuMP.value.(simulation.mechanical_step.prev_θ)
    push!(simulation.steps, SimulationStep(x, y, θ))
end

function simulate!(simulation::Simulation, num_steps=1)
    steps = 1:num_steps
    if num_steps > 1
        steps = ProgressBars.ProgressBar(1:num_steps)
    end

    for _ in steps
        update!(simulation.mechanical_step)
        solve!(simulation.mechanical_step)
        append_step!(simulation)
    end
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
    # CairoMakie.hidedecorations!(ax)
    # CairoMakie.hidespines!(ax)
    CairoMakie.poly!(vertex_coords, faces, color=:transparent, strokewidth=strokewidth, shading=true)

    fig
end

function get_strains(m, prev_x, prev_y, x, y, grid)
    prev_triangles = [[NLExprVector(m, [prev_x[i], prev_y[i]]) for i in T] for T in grid.triangles]
    triangles = [[NLExprVector(m, [x[i], y[i]]) for i in T] for T in grid.triangles]
    reference_triangles = [[[grid.x[i], grid.y[i]] for i in T] for T in grid.triangles]
    prev_strains = strain.(zip(prev_triangles, reference_triangles))
    strains = strain.(zip(triangles, reference_triangles))

    prev_strains, strains
end

function get_symmetrized_strain_rates(prev_strains, strains)
    strain_rates = [F - prev_F for (F, prev_F) in zip(strains, prev_strains)]
    [
        transpose(dot_F) * prev_F + transpose(prev_F) * dot_F
        for (prev_F, dot_F) in zip(prev_strains, strain_rates)
    ]
end

function integrate(f, node_values, m)
    f_symb = Symbol(f)
    θ1, θ2, θ3 = node_values
    expr = :(($f_symb($(θ1)) + $f_symb($(θ2)) + $f_symb($(θ3))) / 3)
    JuMP.add_nonlinear_expression(m, expr)
end

function create_objective!(mechanical_step::MechanicalStep, grid::SimulationGrid, shape_memory_scaling::Number, fps::Number)
    m = mechanical_step.model
    prev_x = mechanical_step.prev_x
    prev_y = mechanical_step.prev_y
    prev_θ = mechanical_step.prev_θ
    x = mechanical_step.x
    y = mechanical_step.y

    # compute strains and symmetrized strain-rates
    prev_strains, strains = get_strains(m, prev_x, prev_y, x, y, grid)
    symmetrized_strain_rates = get_symmetrized_strain_rates(prev_strains, strains)
    JuMP.register(m, :austenite_percentage, 1, austenite_percentage; autodiff=true)
    austenite_percentages = [
        integrate(austenite_percentage, (prev_θ[i1], prev_θ[i2], prev_θ[i3]), m)
        for (i1, i2, i3) in grid.triangles
    ]
    display(JuMP.value.(austenite_percentages))
    # austenite_percentages = Vector{JuMP.NonlinearExpression}(undef, length(grid.triangles))
    # for (i, (i1, i2, i3)) in enumerate(grid.triangles)
    #     austenite_percentages[i] = JuMP.@NLexpression(
    #         m,
    #         (austenite_percentage(prev_θ[i1]) + austenite_percentage(prev_θ[i2]) + austenite_percentage(prev_θ[i3])) / 3
    #     )
    # end

    # objective
    scaling_matrix = [1/shape_memory_scaling 0; 0 1]
    scaled_strains = [F * scaling_matrix for F in strains]
    JuMP.@NLexpression(
        m, elastic_energy,
        0.5 * sum(
            a_perc * neo_hook_F + (1 - a_perc) * neo_hook_F_scaled
            for (F, a_perc, neo_hook_F, neo_hook_F_scaled) in zip(
                strains,
                austenite_percentages,
                neo_hook.(strains),
                neo_hook.(scaled_strains)
            )
        )
    )
    JuMP.@NLexpression(m, dissipation, 0.5 * sum(d for d in norm_sqr.(symmetrized_strain_rates)))
    JuMP.@NLobjective(m, Min, elastic_energy + fps * dissipation)
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
struct SimulationStep
    x::Vector{Float64}
    y::Vector{Float64}
    θ::Vector{Float64}
end

Base.@kwdef struct Simulation
    grid::SimulationGrid
    fps::Number = 30
    initial_temperature::Number = 0.0
    initial_scaling::Number = 1.0
    shape_memory_scaling::Number = 2.0
    temperature_search_radius::Number
    deformation_search_radius::Number
    heat_conductivity::Vector{Number} = [1.0, 1.0]
    heat_transfer_coefficient::Number = 0.5
    entropic_heat_capacity::Number = 10.0
    external_temperature::Number = 0.0
    mechanical_step::JuMP.Model # TODO: create custom struct
    thermal_step::JuMP.Model # TODO: create custom struct
    steps::Vector{SimulationStep}
end

function Simulation(grid::SimulationGrid, deformation_search_radius, temperature_search_radius)
    mechanical_step = create_mechanical_step(grid, deformation_search_radius)
    thermal_step = create_thermal_step()
    steps = Vector{SimulationStep}()
    Simulation(grid=grid, mechanical_step=mechanical_step, thermal_step=thermal_step, steps=steps, deformation_search_radius=deformation_search_radius, temperature_search_radius=temperature_search_radius)
end

function create_mechanical_step(grid::SimulationGrid, search_rad)
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
    JuMP.register(m, :integrate, 3, integrate; autodiff=true)
    austenite_percentages = Vector{JuMP.NonlinearExpression}(undef, length(grid.triangles))
    for (i, (i1, i2, i3)) in enumerate(grid.triangles)
        austenite_percentages[i] = JuMP.@NLexpression(
            m,
            (austenite_percentage(prev_θ[i1]) + austenite_percentage(prev_θ[i2]) + austenite_percentage(prev_θ[i3])) / 3
        )
    end

    m
end

function integrate(f1, f2, f3)
    (f1 + f2 + f3) / 3
end

function austenite_percentage(θ)
    return 1 - 1 / (1 + θ)
end

function strain(triangle, reference_triangle)
    a, b, c = reference_triangle
    gradient_equilateral_to_reference = [b - a c - a]
    a, b, c = triangle
    gradient_equilateral_to_current = [b - a c - a]

    gradient_equilateral_to_current * inv(gradient_equilateral_to_reference)
end

function strain((triangle, reference_triangle))
    strain(triangle, reference_triangle)
end

function create_thermal_step()
    m = JuMP.Model(() -> MadNLP.Optimizer(print_level=MadNLP.WARN, blas_num_threads=4))
    m
end
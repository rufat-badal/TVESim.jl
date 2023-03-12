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
    m = JuMP.Model(() -> MadNLP.Optimizer(print_level=MadNLP.WARN, blas_num_threads=8))

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

mutable struct Simulation
    grid::SimulationGrid
    deformation_search_radius::Number
    temperature_search_radius::Number
    shape_memory_scaling::Number
    initial_temperature::Number
    fps::Number
    heat_transfer_coefficient::Number
    heat_conductivity::Matrix{Number}
    entropic_heat_capacity::Number
    external_temperature::Number
    mechanical_step::MechanicalStep
    thermal_step::ThermalStep
    steps::Vector{SimulationStep}
    x_range::Tuple{Float64,Float64}
    y_range::Tuple{Float64,Float64}
    θ_range::Tuple{Float64,Float64}

    function Simulation(
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
        thermal_step,
        steps
    )
        new(
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
            thermal_step,
            steps,
            get_initial_range(steps[1].x, steps[2].x),
            get_initial_range(steps[1].y, steps[2].y),
            get_initial_range(steps[1].θ, steps[2].θ)
        )
    end
end

get_initial_range(step1, step2) = (
    min(minimum(step1), minimum(step2)),
    max(maximum(step1), maximum(step2))
)

function Simulation(
    grid::SimulationGrid;
    shape_memory_scaling=1.5,
    initial_temperature=0,
    fps=30,
    heat_transfer_coefficient=0.5,
    heat_conductivity=[1 0; 0 1],
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
    create_objective!(
        thermal_step, grid, shape_memory_scaling,
        heat_transfer_coefficient, heat_conductivity, entropic_heat_capacity, external_temperature, fps
    )
    JuMP.optimize!(thermal_step.model)

    x = JuMP.value.(mechanical_step.x)
    y = JuMP.value.(mechanical_step.y)
    θ = JuMP.value.(thermal_step.θ)
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
        thermal_step,
        steps
    )
end

function create_objective!(
    mechanical_step::MechanicalStep,
    grid::SimulationGrid,
    shape_memory_scaling::Number,
    fps
)
    m = mechanical_step.model
    prev_x = mechanical_step.prev_x
    prev_y = mechanical_step.prev_y
    prev_θ = mechanical_step.prev_θ
    x = mechanical_step.x
    y = mechanical_step.y

    # compute strains and symmetrized strain-rates
    prev_strains, strains = get_strains(m, prev_x, prev_y, x, y, grid)
    symmetrized_strain_rates = get_symmetrized_strain_rates(prev_strains, strains)

    # objective
    scaling_matrix = [1/shape_memory_scaling 0; 0 1]
    scaled_strains = [F * scaling_matrix for F in strains]
    JuMP.register(m, :austenite_percentage, 1, austenite_percentage; autodiff=true)
    elastic_energy = add_nonlinear_expression(0.5 * sum(
        a_perc * neo_hook(F) + (1 - a_perc) * neo_hook(F_scaled)
        for (a_perc, F, F_scaled) in zip(
            integral(austenite_percentage, prev_θ, grid, m),
            strains,
            scaled_strains
        )
    ))

    dissipation = add_nonlinear_expression(sum(dissipation_potential(dot_C) for dot_C in symmetrized_strain_rates))
    JuMP.@NLobjective(m, Min, elastic_energy + fps * dissipation)
end

function create_objective!(
    thermal_step::ThermalStep, grid::SimulationGrid,
    shape_memory_scaling,
    heat_transfer_coefficient,
    heat_conductivity,
    entropic_heat_capacity,
    external_temperature,
    fps
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
    strain_rates, symmetrized_strain_rates = get_strain_rates(prev_strains, strains)

    # account for heat transfer on the boundar
    # TODO: use proper approximation of a boundary integral (at least account for the length of each boundary edge)
    boundary_heat_transfer = JuMP.@NLexpression(
        m,
        0.5 * heat_transfer_coefficient * sum(
            (θ[i] - external_temperature)^2 for i in grid.boundary_vertices
        )
    )

    # heat diffusion
    temp_gradients = get_gradients(θ, grid, m)
    heat_diffusion = add_nonlinear_expression(
        0.5 * sum(
            dot(∇θ, heat_conductivity * ∇θ)
            for ∇θ in temp_gradients
        )
    )

    # heat sources and sinks
    scaling_matrix = [1/shape_memory_scaling 0; 0 1]
    JuMP.register(m, :austenite_percentage, 1, austenite_percentage; autodiff=true)
    JuMP.register(m, :identity, 1, identity; autodiff=true)
    adiabatic_terms = [
        dot(
            a_perc_term * (
                gradient_austenite_potential(prev_F) - gradient_martensite_potential(prev_F, scaling_matrix)
            ), dot_F
        )
        for (a_perc_term, prev_F, dot_F) in zip(
            scalar_product(austenite_percentage, prev_θ, θ, grid, m),
            prev_strains,
            strain_rates
        )
    ]

    dissipation_rate(dot_C) = 2 * dissipation_potential(dot_C)
    heat_creation_consumption = add_nonlinear_expression(
        -sum(
            adiab + d_rate * temp
            for (adiab, d_rate, temp) in zip(
                adiabatic_terms,
                dissipation_rate.(symmetrized_strain_rates),
                integral(θ, grid, m),
            )
        )
    )

    # dissipation
    JuMP.register(m, :internal_energy_weight, 1, internal_energy_weight; autodiff=true)
    antider_prev_internal_energies = [
        weight * (
            austenite_potential(prev_F) - martensite_potential(prev_F, scaling_matrix)
        ) + entropic_heat_capacity * entropy_term
        for (weight, prev_F, entropy_term) in zip(
            scalar_product(internal_energy_weight, prev_θ, θ, grid, m),
            prev_strains,
            scalar_product(prev_θ, θ, grid, m)
        )
    ]

    square(x) = x^2
    JuMP.register(m, :square, 1, square; autodiff=true)
    JuMP.register(m, :antider_internal_energy_weight, 1, antider_internal_energy_weight; autodiff=true)
    antider_internal_energies = [
        antider_weight * (
            austenite_potential(F) - martensite_potential(F, scaling_matrix)
        ) + entropic_heat_capacity / 2 * temp_squared
        for (antider_weight, F, temp_squared) in zip(
            integral(antider_internal_energy_weight, θ, grid, m),
            strains,
            integral(square, θ, grid, m)
        )
    ]

    dissipation = add_nonlinear_expression(
        0.5 * sum(
            antider_W - antider_prev_W
            for (antider_W, antider_prev_W) in zip(
                antider_internal_energies,
                antider_prev_internal_energies
            )
        )
    )

    JuMP.@NLobjective(
        m, Min,
        boundary_heat_transfer
        + heat_diffusion
        + heat_creation_consumption
        + fps * dissipation
    )
end

function get_strains(m, prev_x, prev_y, x, y, grid)
    prev_triangles = [[jumpexpression_array(m, [prev_x[i], prev_y[i]]) for i in T] for T in grid.triangles]
    triangles = [[jumpexpression_array(m, [x[i], y[i]]) for i in T] for T in grid.triangles]
    reference_triangles = [[[grid.x[i], grid.y[i]] for i in T] for T in grid.triangles]
    prev_strains = [
        strain(prev_triangle, reference_triangle)
        for (prev_triangle, reference_triangle) in zip(prev_triangles, reference_triangles)
    ]
    strains = [
        strain(triangle, reference_triangle)
        for (triangle, reference_triangle) in zip(triangles, reference_triangles)
    ]

    prev_strains, strains
end

function strain(triangle, reference_triangle)
    a, b, c = reference_triangle
    gradient_equilateral_to_reference = [b - a c - a]
    a, b, c = triangle
    gradient_equilateral_to_current = [b - a c - a]

    gradient_equilateral_to_current * inv(gradient_equilateral_to_reference)
end

function get_symmetrized_strain_rates(prev_strains, strains)
    strain_rates = [F - prev_F for (F, prev_F) in zip(strains, prev_strains)]
    [
        transpose(dot_F) * prev_F + transpose(prev_F) * dot_F
        for (prev_F, dot_F) in zip(prev_strains, strain_rates)
    ]
end

function get_strain_rates(prev_strains, strains)
    strain_rates = [F - prev_F for (F, prev_F) in zip(strains, prev_strains)]
    symmetrized_strain_rates = [
        transpose(dot_F) * prev_F + transpose(prev_F) * dot_F
        for (prev_F, dot_F) in zip(prev_strains, strain_rates)
    ]
    strain_rates, symmetrized_strain_rates
end

function get_gradients(f, grid, model)
    node_values = [jumpexpression_array(model, [f[i1], f[i2], f[i3]]) for (i1, i2, i3) in grid.triangles]
    reference_triangles = [[[grid.x[i], grid.y[i]] for i in T] for T in grid.triangles]
    [gradient(f, reference_triangle) for (f, reference_triangle) in zip(node_values, reference_triangles)]
end

function gradient(f, reference_triangle)
    fa, fb, fc = f
    a, b, c = reference_triangle

    transpose([fb - fa fc - fa] * inv([b - a c - a]))
end

function integral(f::Function, x, grid, model::JuMP.Model)
    # It is assumed that f was already registered by the caller
    [
        integral(f, (x[i1], x[i2], x[i3]), model) * a_fac
        for ((i1, i2, i3), a_fac) in zip(
            grid.triangles,
            grid.area_factors
        )
    ]
end

function integral(f::Function, x, model::JuMP.Model)
    f_symb = Symbol(f)
    x1, x2, x3 = x
    expr = :(($f_symb($(x1)) + $f_symb($(x2)) + $f_symb($(x3))) / 6)
    JuMPExpression(model, expr)
end

integral(x, grid, model) = integral(identity, x, grid, model)

function scalar_product(f::Function, x, g::Function, y, grid, model::JuMP.Model)
    # It is assumed that f and g were already registered by the caller
    [
        scalar_product(f, (x[i1], x[i2], x[i3]), g, (y[i1], y[i2], y[i3]), model) * a_fac
        for ((i1, i2, i3), a_fac) in zip(
            grid.triangles,
            grid.area_factors
        )
    ]
end

function scalar_product(f::Function, x, g::Function, y, model::JuMP.Model)
    f_symb = Symbol(f)
    x1, x2, x3 = x
    g_symb = Symbol(g)
    y1, y2, y3 = y

    expr = :(
        (
            $f_symb($(x1)) * $g_symb($(y1))
            + $f_symb($(x2)) * $g_symb($(y2))
            + $f_symb($(x3)) * $g_symb($(y3))
        ) / 6
    )
    JuMPExpression(model, expr)
end

scalar_product(x, y, grid, model::JuMP.Model) = scalar_product(identity, x, identity, y, grid, model)
scalar_product(f::Function, x, y, grid, model::JuMP.Model) = scalar_product(f, x, identity, y, grid, model)
scalar_product(x, g::Function, y, grid, model::JuMP.Model) = scalar_product(identity, x, g, y, grid, model)

function update_mechanical_step!(simulation::Simulation)
    mechanical_step = simulation.mechanical_step
    thermal_step = simulation.thermal_step
    JuMP.set_value.(mechanical_step.prev_θ, JuMP.value.(thermal_step.θ))
    JuMP.set_value.(mechanical_step.prev_x, JuMP.value.(mechanical_step.x))
    JuMP.set_value.(mechanical_step.prev_y, JuMP.value.(mechanical_step.y))
end

function update_thermal_step!(simulation::Simulation)
    mechanical_step = simulation.mechanical_step
    thermal_step = simulation.thermal_step
    JuMP.set_value.(thermal_step.prev_x, JuMP.value.(mechanical_step.prev_x))
    JuMP.set_value.(thermal_step.prev_y, JuMP.value.(mechanical_step.prev_y))
    JuMP.set_value.(thermal_step.x, JuMP.value.(mechanical_step.x))
    JuMP.set_value.(thermal_step.y, JuMP.value.(mechanical_step.y))
    JuMP.set_value.(thermal_step.prev_θ, JuMP.value.(thermal_step.θ))
end

function solve!(mechanical_step::MechanicalStep)
    JuMP.optimize!(mechanical_step.model)
end

function solve!(thermal_step::ThermalStep)
    JuMP.optimize!(thermal_step.model)
end

function update_ranges!(simulation)
    simulation.x_range = get_new_range(simulation.x_range, simulation.steps[end].x)
    simulation.y_range = get_new_range(simulation.y_range, simulation.steps[end].y)
    simulation.θ_range = get_new_range(simulation.θ_range, simulation.steps[end].θ)
end

get_new_range(old_range, new_step) = (
    min(old_range[1], minimum(new_step)),
    max(old_range[2], maximum(new_step))
)

function append_step!(simulation::Simulation)
    x = JuMP.value.(simulation.mechanical_step.x)
    y = JuMP.value.(simulation.mechanical_step.y)
    θ = JuMP.value.(simulation.thermal_step.θ)
    push!(simulation.steps, SimulationStep(x, y, θ))
    update_ranges!(simulation)
end

function simulate!(simulation::Simulation, num_steps=1)
    steps = 1:num_steps
    if num_steps > 1
        steps = ProgressBars.ProgressBar(1:num_steps)
    end

    for _ in steps
        update_mechanical_step!(simulation)
        solve!(simulation.mechanical_step)
        update_thermal_step!(simulation)
        solve!(simulation.thermal_step)
        append_step!(simulation)
    end
end

function plot(simulation::Simulation, i::Int; show_edges=false)
    num_horizontal_pixels = 2500
    strokewidth = 1 / 1000 * num_horizontal_pixels

    step = simulation.steps[i]
    x, y = step.x, step.y
    vertex_coords = [x y]
    triangles = simulation.grid.triangles
    # TODO: generate this in the grid object
    faces = Matrix{Int}(undef, length(triangles), 3)
    for (i, T) in enumerate(triangles)
        for (j, v) in enumerate(T)
            faces[i, j] = v
        end
    end
    vertices = [x y]

    min_x, max_x = simulation.x_range
    min_y, max_y = simulation.y_range
    width = max_x - min_x
    height = max_y - min_y
    aspect = width / height
    plot_height = num_horizontal_pixels / aspect
    fig = CairoMakie.Figure(resolution=(num_horizontal_pixels, plot_height))
    padding = 0
    if show_edges
        length_pixel = width / num_horizontal_pixels
        padding = strokewidth / 2 * length_pixel
    end
    ax = CairoMakie.Axis(
        fig[1, 1],
        limits=(
            min_x - padding, max_x + padding,
            min_y - padding, max_y + padding),
        aspect=aspect)
    CairoMakie.hidedecorations!(ax)
    CairoMakie.hidespines!(ax)

    if show_edges
        CairoMakie.poly!(
            vertices,
            faces,
            color=step.θ,
            colormap=:plasma,
            colorrange=simulation.θ_range,
            strokewidth=strokewidth,
            shading=true
        )
    else
        CairoMakie.mesh!(
            vertices,
            faces,
            color=step.θ,
            colormap=:plasma,
            colorrange=simulation.θ_range
        )
    end

    fig
end

function save(simulation::Simulation, folder; show_edges=false)
    for i in 1:length(simulation.steps)
        CairoMakie.save("$folder/step_$i.png", TVESim.plot(simulation, i, show_edges=show_edges))
    end
end

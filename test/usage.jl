using TVESim
using CairoMakie
using Dates

function setup_experiment_directory(keyword)
    t = Dates.format(now(), "yy-mm-dd-H-M-S")
    mkpath("results/experiments/$(keyword)_$t/")
end

function square_experiment()
    experiments_directory = setup_experiment_directory("square")
    width = 10
    height = 10
    initial_temperature = 0.0
    isinternal(x, y) = -(width + 1) / 2 <= x <= (width + 1) / 2 && -(height + 1) / 2 <= y <= (height + 1) / 2
    isdirichlet(x, y) = x <= -(width - 1) / 2
    grid = SimulationGrid(width + 2, height + 2, initial_temperature, TVESim.isosceles_right_triangulation, isinternal, isdirichlet)
    simulation = Simulation(grid, fps=0.75, heat_transfer_coefficient=6)
    num_steps = 40
    simulate!(simulation, num_steps - 2)
    for i in 1:num_steps
        save("$experiments_directory/step_$i.png", TVESim.plot(simulation, i, show_edges=true))
    end
    save("$experiments_directory/grid.png", TVESim.plot(grid, simulation.θ_range, show_edges=true))
    return
end

function cooling_square_experiment()
    experiments_directory = setup_experiment_directory("cooling_square")
    width = 10
    height = 10
    initial_temperature = 10
    isinternal(x, y) = -(width + 1) / 2 <= x <= (width + 1) / 2 && -(height + 1) / 2 <= y <= (height + 1) / 2
    isdirichlet(x, y) = x <= -(width - 1) / 2
    grid = SimulationGrid(width + 2, height + 2, initial_temperature, TVESim.isosceles_right_triangulation, isinternal, isdirichlet)
    simulation = Simulation(grid, fps=1.5, heat_transfer_coefficient=4)
    num_steps = 50
    simulate!(simulation, num_steps - 2)
    for i in 1:num_steps
        save("$experiments_directory/step_$i.png", TVESim.plot(simulation, i, show_edges=true))
    end
    save("$experiments_directory/grid.png", TVESim.plot(grid, simulation.θ_range, show_edges=true))
    return
end

function circle_experiment()
    experiments_directory = setup_experiment_directory("circle")
    initial_temperature = 0.0
    radius = 5
    num_boundary_points = 70
    dirichlet_arc_angle = 45
    isdirichlet(x, y) = x <= -radius * cos(dirichlet_arc_angle / 360 * 2pi)
    grid = SimulationGrid(TVESim.circle_boundary_points(radius, num_boundary_points), initial_temperature, isdirichlet)
    simulation = Simulation(grid, fps=1.5, heat_transfer_coefficient=6)
    num_steps = 20
    simulate!(simulation, num_steps - 2)
    for i in 1:num_steps
        save("$experiments_directory/step_$i.png", TVESim.plot(simulation, i, show_edges=true))
    end
    save("$experiments_directory/grid.png", TVESim.plot(grid, simulation.θ_range, show_edges=true))
    return
end

circle_experiment()
square_experiment()
cooling_square_experiment()
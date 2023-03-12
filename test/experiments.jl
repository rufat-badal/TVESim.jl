using TVESim
using Dates

function setup_experiment_folder(keyword)
    t = Dates.format(now(), "dd-mm-yyyy-H-M-S")
    mkpath("results/experiments/$(keyword)_$t/")
end

function square_experiment()
    experiments_folder = setup_experiment_folder("square")
    width = 10
    height = 10
    initial_temperature = 0.0
    isinternal(x, y) = -(width + 1) / 2 <= x <= (width + 1) / 2 && -(height + 1) / 2 <= y <= (height + 1) / 2
    isdirichlet(x, y) = x <= -(width - 1) / 2
    grid = SimulationGrid(width + 2, height + 2, initial_temperature, TVESim.isosceles_right_triangulation, isinternal, isdirichlet)
    simulation = Simulation(grid, fps=0.75, heat_transfer_coefficient=6)
    num_steps = 40
    simulate!(simulation, num_steps - 2)
    save(simulation, experiments_folder, show_edges=true)
    save(grid, simulation.θ_range, experiments_folder, show_edges=true)
    return
end

function cooling_square_experiment()
    experiments_folder = setup_experiment_folder("cooling_square")
    width = 10
    height = 10
    initial_temperature = 10
    isinternal(x, y) = -(width + 1) / 2 <= x <= (width + 1) / 2 && -(height + 1) / 2 <= y <= (height + 1) / 2
    isdirichlet(x, y) = x <= -(width - 1) / 2
    grid = SimulationGrid(width + 2, height + 2, initial_temperature, TVESim.isosceles_right_triangulation, isinternal, isdirichlet)
    simulation = Simulation(grid, fps=1.3, heat_transfer_coefficient=4)
    num_steps = 80
    simulate!(simulation, num_steps - 2)
    save(simulation, experiments_folder, show_edges=true)
    save(grid, simulation.θ_range, experiments_folder, show_edges=true)
    return
end

function circle_experiment()
    experiments_folder = setup_experiment_folder("circle")
    initial_temperature = 0.0
    radius = 5
    num_boundary_points = 70
    dirichlet_arc_angle = 45
    isdirichlet(x, y) = x <= -radius * cos(dirichlet_arc_angle / 360 * 2pi)
    grid = SimulationGrid(TVESim.circle_boundary_points(radius, num_boundary_points), initial_temperature, isdirichlet)
    simulation = Simulation(grid, fps=1.5, heat_transfer_coefficient=6)
    num_steps = 20
    simulate!(simulation, num_steps - 2)
    save(simulation, experiments_folder, show_edges=true)
    save(grid, simulation.θ_range, experiments_folder, show_edges=true)
    return
end

square_experiment()
# cooling_square_experiment()
# circle_experiment()
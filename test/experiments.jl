using TVESim
using Dates

function setup_experiment_folder(keyword)
    t = Dates.format(now(), "dd-mm-yyyy-H-M-S")
    mkpath("results/experiments/$(keyword)_$t/")
end

function square_movie_experiment()
    experiments_folder = setup_experiment_folder("square")
    width = 10
    height = 10
    initial_temperature = 0.0
    isinternal(x, y) = -(width + 1) / 2 <= x <= (width + 1) / 2 && -(height + 1) / 2 <= y <= (height + 1) / 2
    isdirichlet(x, y) = x <= -(width - 1) / 2
    grid = SimulationGrid(width + 2, height + 2, initial_temperature, TVESim.isosceles_right_triangulation, isinternal, isdirichlet)
    simulation = Simulation(grid, fps=1, heat_transfer_coefficient=6)
    num_steps = 40
    simulate!(simulation, num_steps - 2)
    save(grid, simulation.θ_range, experiments_folder, show_edges=true)
    save(simulation, experiments_folder, show_edges=true)
    return
end

function square_movie_experiment()
    experiments_folder = setup_experiment_folder("square")
    width = 10
    height = 10
    initial_temperature = 0.0
    isinternal(x, y) = -(width + 1) / 2 <= x <= (width + 1) / 2 && -(height + 1) / 2 <= y <= (height + 1) / 2
    isdirichlet(x, y) = x <= -(width - 1) / 2
    grid = SimulationGrid(width + 2, height + 2, initial_temperature, TVESim.isosceles_right_triangulation, isinternal, isdirichlet)
    simulation = Simulation(grid, fps=3, heat_transfer_coefficient=6)
    num_steps = 90
    simulate!(simulation, num_steps - 2)
    save(grid, simulation.θ_range, experiments_folder, show_edges=false)
    save(simulation, experiments_folder, show_edges=true, movie=true)
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
    simulation = Simulation(grid, fps=1, heat_transfer_coefficient=4)
    num_steps = 80
    simulate!(simulation, num_steps - 2)
    save(grid, simulation.θ_range, experiments_folder, show_edges=true)
    save(simulation, experiments_folder, show_edges=true)
    return
end

function circle_experiment()
    experiments_folder = setup_experiment_folder("circle")
    initial_temperature = 0.0
    radius = 5
    num_boundary_points = 30
    dirichlet_arc_angle = 45
    isdirichlet(x, y) = x <= -radius * cos(dirichlet_arc_angle / 360 * 2pi)
    grid = SimulationGrid(TVESim.circle_boundary_points(radius, num_boundary_points), initial_temperature, isdirichlet)
    simulation = Simulation(grid, fps=1.5, heat_transfer_coefficient=6)
    num_steps = 20
    simulate!(simulation, num_steps - 2)
    save(grid, simulation.θ_range, experiments_folder, show_edges=true)
    save(simulation, experiments_folder, show_edges=true)
    return
end

function cooling_circle_experiment()
    experiments_folder = setup_experiment_folder("cooling_circle")
    initial_temperature = 10
    radius = 5
    num_boundary_points = 30
    dirichlet_arc_angle = 45
    isdirichlet(x, y) = x <= -radius * cos(dirichlet_arc_angle / 360 * 2pi)
    grid = SimulationGrid(TVESim.circle_boundary_points(radius, num_boundary_points), initial_temperature, isdirichlet)
    simulation = Simulation(grid, fps=1, heat_transfer_coefficient=4)
    num_steps = 80
    simulate!(simulation, num_steps - 2)
    save(grid, simulation.θ_range, experiments_folder, show_edges=true)
    save(simulation, experiments_folder, show_edges=true)
    return
end

function rotating_circle_experiment()
    rotation_time = 1
    experiments_folder = setup_experiment_folder("rotating_circle")
    initial_temperature = 5
    radius = 5
    num_boundary_points = 30
    isdirichlet(x, y) = true
    function boundary_rotation(x, y, t)
        angle = 2 * π * t / rotation_time
        c = cos(angle)
        s = sin(angle)
        c * x - s * y, s * x + c * y
    end
    grid = SimulationGrid(TVESim.circle_boundary_points(radius, num_boundary_points), initial_temperature, isdirichlet)
    simulation = Simulation(
        grid,
        dirichlet_func=boundary_rotation,
        fps=30,
        heat_transfer_coefficient=0
    )
    num_steps = 30
    simulate!(simulation, num_steps - 2)
    save(grid, simulation.θ_range, experiments_folder, show_edges=true)
    save(simulation, experiments_folder, show_edges=true)
    return simulation.θ_range
end

# square_experiment()
# cooling_square_experiment()
# circle_experiment()
# cooling_circle_experiment()
# square_movie_experiment()
rotating_circle_experiment()
using TVESim
using CairoMakie

radius = 20
initial_temperature = 0.0
num_boundary_points = 25
dirichlet_arc_angle = 45
isdirichlet(x, y) = x - radius <= -radius * cos(dirichlet_arc_angle / 360 * 2pi)
grid = SimulationGrid(TVESim.circle_boundary_points(radius, num_boundary_points), initial_temperature, isdirichlet)
display(TVESim.plot(grid, (0.0, 1.0), show_edges=true))
sim = Simulation(grid, 1.0, 5.0)
sim.mechanical_step;
using TVESim
using CairoMakie

initial_temperature = 0.0
radius = 20
num_boundary_points = 50
dirichlet_arc_angle = 45
isdirichlet(x, y) = x - radius <= -radius * cos(dirichlet_arc_angle / 360 * 2pi)
# isdirichlet(x, y) = false
grid = SimulationGrid(TVESim.circle_boundary_points(radius, num_boundary_points), initial_temperature, isdirichlet)
# save("results/grid.png", TVESim.plot(grid, (0, 1), show_edges=true))
simulation = Simulation(grid, fps=10, heat_transfer_coefficient=0)
# width = 10
# height = 10
# isinternal(x, y) = 1 <= x <= width + 1 && 1 <= y <= height + 1
# isdirichlet(x, y) = x <= 1.5
# grid = SimulationGrid(width + 2, height + 2, initial_temperature, TVESim.isosceles_right_triangulation, isinternal, isdirichlet)
num_steps=40
for i in 3:num_steps
    simulate!(simulation)
end
for i in 1:num_steps
    save("results/step_$i.png", TVESim.plot(simulation, i, show_edges=true))
end
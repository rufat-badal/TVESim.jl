using TVESim
using CairoMakie

initial_temperature = 0.0
radius = 20
num_boundary_points = 50
dirichlet_arc_angle = 45
isdirichlet(x, y) = x - radius <= -radius * cos(dirichlet_arc_angle / 360 * 2pi)
# isdirichlet(x, y) = false
grid = SimulationGrid(TVESim.circle_boundary_points(radius, num_boundary_points), initial_temperature, isdirichlet)
save("results/grid.png", TVESim.plot(grid, (0, 1), show_edges=true))
simulation = Simulation(grid, fps=5, heat_transfer_coefficient=0)
display(simulation.x_range)
display(simulation.θ_range)
# width = 10
# height = 10
# isinternal(x, y) = 1 <= x <= width + 1 && 1 <= y <= height + 1
# isdirichlet(x, y) = x <= 1.5
# grid = SimulationGrid(width + 2, height + 2, initial_temperature, TVESim.isosceles_right_triangulation, isinternal, isdirichlet)
save("results/step_1.png", TVESim.plot(simulation.steps[1], simulation.grid.triangles))
save("results/step_2.png", TVESim.plot(simulation.steps[2], simulation.grid.triangles))
num_steps=10
for i in 3:num_steps
    simulate!(simulation)
    save("results/step_$i.png", TVESim.plot(simulation.steps[end], simulation.grid.triangles))
end
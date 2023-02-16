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
    temp_search_radius::Number = 1.0
    heat_conductivity::Vector{Number} = [1.0, 1.0]
    heat_transfer_coefficient::Number = 0.5
    entropic_heat_capacity::Number = 10.0
    external_temperature::Number = 0.0
    mechanical_step::JuMP.Model
    thermal_step::JuMP.Model
    steps::Vector{SimulationStep}
end

function Simulation(grid::SimulationGrid)
    mechanical_step = JuMP.Model(() -> MadNLP.Optimizer(print_level=MadNLP.WARN, blas_num_threads=8))
    thermal_step = JuMP.Model(() -> MadNLP.Optimizer(print_level=MadNLP.WARN, blas_num_threads=4))
    steps = Vector{SimulationStep}()
    Simulation(grid=grid, mechanical_step=mechanical_step, thermal_step=thermal_step, steps=steps)
end
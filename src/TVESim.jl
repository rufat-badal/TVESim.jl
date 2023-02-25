module TVESim

export SimulationGrid
export Simulation
export simulate!

import JuMP
import JuMP.value
import MadNLP
import Base
import Base: +, -, *, inv
import LinearAlgebra: tr, det, transpose
import CairoMakie
import Triangulate
import ProgressBars

include("linalg_for_nlexpr.jl")
include("grid.jl")
include("simulation.jl")

end

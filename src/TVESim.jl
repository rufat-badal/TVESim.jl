module TVESim

export SimulationGrid
export Simulation

import JuMP
import JuMP.value
import MadNLP
import Base
import Base: +, -, *, inv
import LinearAlgebra: tr, det, transpose
import CairoMakie
import Triangulate

include("linalg_for_nlexpr.jl")
include("grid.jl")
include("simulation.jl")

end

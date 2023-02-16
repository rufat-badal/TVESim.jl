module TVESim

export SimulationGrid
export plot

import JuMP
import JuMP.value
import Base: +, -, *
import LinearAlgebra: tr, det, transpose
import CairoMakie
import Triangulate

include("linalg_for_nlexpr.jl")
include("grid.jl")

end

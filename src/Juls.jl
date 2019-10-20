module Juls

export RunJuls, Parameter

#using NetCDF
#using Dates, FileIO, Printf, Parameters
using Parameters, Printf, Dates

include("DefaultParameters.jl")
include("Grid.jl")
include("Constants.jl")
include("Forcing.jl")
include("ModelSetup.jl")
include("InitialConditions.jl")
include("Preallocate.jl")

include("TimeIntegration.jl")
include("GhostPoints.jl")
include("rhs.jl")
include("Gradients.jl")
include("Interpolations.jl")
include("PVadvection.jl")
include("Continuity.jl")
include("Bottomdrag.jl")
include("Diffusion.jl")
#include("TracerAdvection.jl")

include("Feedback.jl")
#include("Output.jl")
include("RunJuls.jl")

end

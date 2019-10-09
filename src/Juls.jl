module Juls

export RunJuls, Parameter

using Dates, NetCDF, FileIO, Printf, Parameters

include("DefaultParameters.jl")
include("Grid.jl")
include("Constants.jl")
include("InitialConditions.jl")
include("Preallocate.jl")

include("TimeIntegration.jl")
include("Forcing.jl")
include("GhostPoints.jl")
include("rhs.jl")
include("Gradients.jl")
include("Interpolations.jl")
include("PVadvection.jl")
include("Continuity.jl")
include("Bottomdrag.jl")
include("Diffusion.jl")
include("TracerAdvection.jl")

include("Feedback.jl")
include("Output.jl")
include("RunJuls.jl")

end

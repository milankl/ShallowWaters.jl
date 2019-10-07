module Juls

export RunJuls, Parameter

using Dates, NetCDF, FileIO, Printf, Parameters

include("DefaultParameters.jl")
include("Grid.jl")
include("RunJuls.jl")
include("Constants.jl")
include("InitialConditions.jl")
include("Preallocate.jl")

# OPERATORS and everything that is needed for the RHS
include("Gradients.jl")
include("Interpolations.jl")

include("PV_adv.jl")
include("bottom_friction.jl")
include("diffusion.jl")
include("TracerAdvection.jl")
include("Forcing.jl")
include("rhs.jl")
include("continuity.jl")
include("time_integration.jl")
include("ghost_points.jl")
include("initial_conditions.jl")
include("preallocate.jl")

# OUTPUT AND FEEDBACK
include("feedback.jl")
include("output.jl")

include("RunJuls.jl")

end

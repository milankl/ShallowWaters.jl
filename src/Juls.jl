module Juls

export RunJuls

using Dates, NetCDF, FileIO, Printf, Parameters

# GRID
include("DefaultParameters.jl")
include("Grid.jl")
include("RunJuls.jl")
include("Constants.jl")


# OPERATORS and everything that is needed for the RHS
include("Gradients.jl")
include("interpolations.jl")
include("PV_adv.jl")
include("bottom_friction.jl")
include("diffusion.jl")
include("tracer_adv.jl")
include("coriolis.jl")
include("forcing.jl")
include("bottom_topography.jl")
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

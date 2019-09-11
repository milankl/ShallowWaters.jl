module Juls

export RunJuls

using Dates, NetCDF, FileIO, Printf

# GRID
include("grid.jl")
include("constants.jl")
#include("domain_decomposition.jl")

# OPERATORS and everything that is needed for the RHS
include("gradients.jl")
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
include("src/feedback.jl")
include("src/output.jl")

end

module ShallowWaters

export RunModel, Parameter, ∂x, ∂y, Ix, Iy, ∇²

using NetCDF, Parameters, Printf, Dates, Interpolations

include("default_parameters.jl")
include("grids.jl")
include("constants.jl")
include("forcing.jl")
include("model_setup.jl")
include("initial_conditions.jl")
include("preallocate.jl")

include("time_integration.jl")
include("ghost_points.jl")
include("rhs.jl")
include("gradients.jl")
include("interpolations.jl")
include("PVadvection.jl")
include("continuity.jl")
include("bottom_drag.jl")
include("diffusion.jl")
include("tracer_advection.jl")

include("feedback.jl")
include("output.jl")
include("run_model.jl")

end

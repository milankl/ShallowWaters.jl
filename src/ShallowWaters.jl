module ShallowWaters

export RunModel, Parameter, ∂x, ∂y, Ix, Iy, ∇²

using NetCDF, Parameters, Printf, Dates, Interpolations

include("default_parametersTEST.jl")
include("gridTEST.jl")
include("constantsTEST.jl")
include("forcingTEST.jl")
include("model_setupTEST.jl")
include("initial_conditionsTEST.jl")
include("preallocateTEST.jl")

include("time_integrationTEST.jl")
include("ghost_pointsTEST.jl")
include("rhsTEST.jl")
include("gradientsTEST.jl")
include("interpolationsTEST.jl")
include("advectionTEST.jl")
include("continuityTEST.jl")
include("bottom_dragTEST.jl")
include("diffusionTEST.jl")
include("tracer_advectionTEST.jl")

include("feedbackTEST.jl")
include("outputTEST.jl")
include("run_modelTEST.jl")

end

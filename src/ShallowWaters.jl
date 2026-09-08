module ShallowWaters

    export run_model, model_setup, Parameter, ∂x, ∂y, Ix, Iy, ∇²

    using NetCDF, Parameters, Printf, Dates, Interpolations

    include("default_parameters.jl")
    include("grid.jl")
    include("constants.jl")
    include("forcing.jl")
    include("preallocate.jl")
    include("model_setup.jl")
    include("initial_conditions.jl")

    include("time_integration.jl")
    include("ghost_points.jl")
    include("rhs.jl")
    include("gradients.jl")
    include("interpolations.jl")
    include("PV_advection.jl")
    include("continuity.jl")
    include("bottom_drag.jl")
    include("diffusion.jl")
    include("tracer_advection.jl")
    include("zanna_bolton_forcing.jl")

    include("feedback.jl")
    include("output.jl")
    include("run_model.jl")

end

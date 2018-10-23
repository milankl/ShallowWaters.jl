#= This script includes all modules and functions, sets up the grid parameters
etc, and runs the model =#

using Dates
using Printf
using NetCDF
using MPI
#using SigmoidNumbers

# PARAMETERS, GRID, CONSTANTS and DOMAIN DECOMPOSITION
include("parameters.jl")
include("src/grid.jl")
include("src/constants.jl")
include("src/domain_decomposition.jl")

# OPERATORS and everything that is needed for the RHS
include("src/gradients.jl")
include("src/interpolations.jl")
include("src/arakawahsu.jl")
include("src/coriolis.jl")
include("src/forcing.jl")
include("src/bottom_topography.jl")
include("src/viscosity.jl")
include("src/rhs.jl")
include("src/time_integration.jl")
include("src/ghost_points.jl")

# OUTPUT AND FEEDBACK
include("src/feedback.jl")
include("src/output.jl")

# INITIALISE
include("src/initial_conditions.jl")
include("src/preallocate.jl")
u,v,η = initial_conditions()

# RUN
u,v,η = time_integration(u,v,η)

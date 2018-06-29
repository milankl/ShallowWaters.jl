#= This script includes all modules and functions, sets up the grid parameters
etc, and runs the model =#

#using JLD
#using netCDF4 at some point in the future for output
using SigmoidNumbers

# PARAMETERS, GRID and CONSTANTS
include("parameters.jl")
include("src/grid.jl")
include("src/constants.jl")

# OPERATORS and everything that is needed for the RHS
include("src/gradients.jl")
include("src/interpolations.jl")
include("src/laplace.jl")
include("src/boundary_conditions.jl")
include("src/coriolis.jl")
include("src/forcing.jl")
include("src/viscosity.jl")
include("src/rhs.jl")
include("src/time_integration.jl")

#TODO
#include("src/nc_output.jl")
#include("src/jld_output.jl")
#include("src/feeback.jl")

# INITILIASE
include("src/initial_conditions.jl")
include("src/preallocate.jl")
u,v,η = initial_conditions()
u,v,η = time_integration(u,v,η)

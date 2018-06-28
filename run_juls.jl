#= This script imports all modules and functions, sets up the grid parameters
etc, and runs the model =#

using JLD
#using netCDF4 at some point in the future for output
#using SigmoidNumbers

# PARAMETERS, GRID and CONSTANTS
import("parameters.jl")
import("src/grid.jl")
import("src/constants.jl")

# OPERATORS and everything that is needed for the RHS
import("src/gradients.jl")
import("src/interpolations.jl")
import("src/laplace.jl")
import("src/boundary_conditions.jl")
import("src/forcing.jl")
import("src/viscosity.jl")
import("src/rhs.jl")
import("src/time_integration.jl")

#TODO
#import("src/nc_output.jl")
#import("src/jld_output.jl")
#import("src/feeback.jl")

# INITILIASE
import("src/initial_conditions.jl")
import("src/preallocate.jl")

u,v,η = initial_conditions()
v_u,dηdx,dLu = preallocate_u_vars(u)
u_v,dηdy,dLv = preallocate_v_vars(v)
dudx,dvdy = preallocate_T_variables(η)

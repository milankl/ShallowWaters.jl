#= This script imports all modules and functions, sets up the grid parameters
etc, and runs the model =#

using JLD
#using netCDF4 at some point in the future for output

import("operators/gradients.jl")
import("operators/interpolations.jl")
import("operators/laplace.jl")

#TODO
#import("integration/rhs.jl")
#import("integration/time_stepping.jl")

#import("output/nc_output.jl")
#import("output/jld_output.jl")
#import("output/feeback.jl")

#TODO
#import("parameters.jl")
# somehow obtain the initial conditions
# preallocate a lot of memory
# start

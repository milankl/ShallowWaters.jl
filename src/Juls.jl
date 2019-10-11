module Juls

export RunJuls, Parameter

#using NetCDF
#using Dates, FileIO, Printf, Parameters
using Parameters

include("DefaultParameters.jl")
include("Grid.jl")
include("Constants.jl")
include("InitialConditions.jl")
include("Preallocate.jl")

include("TimeIntegration.jl")
include("Forcing.jl")
include("GhostPoints.jl")
include("rhs.jl")
include("Gradients.jl")
include("Interpolations.jl")
include("PVadvection.jl")
include("Continuity.jl")
include("Bottomdrag.jl")
include("Diffusion.jl")
#include("TracerAdvection.jl")

#include("Feedback.jl")
#include("Output.jl")
include("RunJuls.jl")

end

using BenchmarkTools
T = Float32
P = Parameter()
G = Grid{T}(P)
C = Constants{T}(P,G)
Prog = initial_conditions(T,P,C,G)
Diag = preallocate(T,G)
Forc = Forcing{T}(P,G)

@unpack u,v,η = Prog
@btime rhs!($u,$v,$η,$P,$C,$G,$Diag,$Forc)

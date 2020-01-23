"""

    u,v,η,sst = RunModel()

runs ShallowWaters with default parameters as defined in src/DefaultParameters.jl

# Examples
```jldoc
julia> u,v,η,sst = RunModel(Float64,nx=200,output=true)
```
"""
function RunModel(::Type{T}=Float32;     # number format
    kwargs...                           # all additional parameters
    ) where {T<:AbstractFloat}

    P = Parameter(T=T;kwargs...)
    return RunModel(T,P)
end

function RunModel(P::Parameter)
    @unpack T = P
    return RunModel(T,P)
end

function RunModel(::Type{T},P::Parameter) where {T<:AbstractFloat}

    @unpack Tprog = P

    G = Grid{T,Tprog}(P)
    C = Constants{T,Tprog}(P,G)
    F = Forcing{T}(P,G)
    S = ModelSetup{T,Tprog}(P,G,C,F)

    Prog = initial_conditions(Tprog,S)
    Diag = preallocate(T,Tprog,G)

    Prog = time_integration(Prog,Diag,S)
    return Prog
end

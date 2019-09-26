"""

    u,v,η,sst = RunJuls()

runs Juls with default parameters as defined in src/DefaultParameters.jl

# Examples
```jldoc
julia> u,v,η,sst = RunJuls(Float64,nx=200,output=true)
```
"""
function RunJuls(::Type{T}=Float32;     # number format
    kwargs...                           # all additional parameters
    ) where {T<:AbstractFloat}

    P = Parameter(T=T;kwargs...)
    G = Grid{T}(P)
    C = Constants{T}(P,G)
    Prog = InitialConditions(T,P,G)
    # Diag = PreallocateDiagnosticVars()
    # S = State(P,G,C,Prog,Diag)
    # time_integration(S)

    return P,G,C,Prog
end

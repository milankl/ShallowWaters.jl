"""

    u,v,η,sst = run_model()

runs ShallowWaters with default parameters as defined in src/DefaultParameters.jl

# Examples
```jldoc
julia> u,v,η,sst = run_model(Float64,nx=200,output=true)
```
"""
function run_model(::Type{T}=Float32;     # number format
    kwargs...                             # all additional parameters
    ) where {T<:AbstractFloat}

    P = Parameter(T=T;kwargs...)
    return run_model(T,P)
end

function run_model(P::Parameter)
    @unpack T = P
    return run_model(T,P)
end

function run_model(::Type{T},P::Parameter) where {T<:AbstractFloat}

    @unpack Tprog = P

    G = Grid{T,Tprog}(P)
    C = Constants{T,Tprog}(P,G)
    F = Forcing{T}(P,G)

    Prog = initial_conditions(Tprog,G,P,C)
    Diag = preallocate(T,Tprog,G)

    # one structure with everything already inside 
    S = ModelSetup{T,Tprog}(P,G,C,F,Prog,Diag,0)
    Prog = time_integration(S)

    return Prog

end

function run_setup(::Type{T}=Float32;     # number format
    kwargs...                             # all additional parameters
    ) where {T<:AbstractFloat}

    P = ShallowWaters.Parameter(T=T;kwargs...)
    return run_setup(T,P)
end

function run_setup(P::ShallowWaters.Parameter)
    @unpack T = P
    return run_setup(T,P)
end

function run_setup(::Type{T},P::ShallowWaters.Parameter) where {T<:AbstractFloat}

    @unpack Tprog = P

    G = ShallowWaters.Grid{T,Tprog}(P)
    C = ShallowWaters.Constants{T,Tprog}(P,G)
    F = ShallowWaters.Forcing{T}(P,G)

    Prog = ShallowWaters.initial_conditions(Tprog,G,P,C)
    Diag = ShallowWaters.preallocate(T,Tprog,G)

    # one structure with everything inside 
    S = ShallowWaters.ModelSetup{T,Tprog}(P,G,C,F,Prog,Diag,0)

    return S

end
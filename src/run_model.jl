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
    # S = ModelSetup{T,Tprog}(P,G,C,F)

    Prog = initial_conditions(Tprog,G,P,C)
    Diag = preallocate(T,Tprog,G)

    # one structure with everything already inside 
    S = ModelSetup{T,Tprog}(P,G,C,F,Prog,Diag,0)
    P = time_integration(S)

    return S, P

end


#### Original ##################################################################################

# In order to run the following you need to change lines in model_setup.jl

function mk_run_model(::Type{T}=Float32;     # number format
    kwargs...                            # all additional parameters
    ) where {T<:AbstractFloat}

    P = Parameter(T=T;kwargs...)
    return mk_run_model(T,P)
end

function mk_run_model(P::Parameter)
    @unpack T = P
    return mk_run_model(T,P)
end

function mk_run_model(::Type{T},P::Parameter) where {T<:AbstractFloat}

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
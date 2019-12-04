struct ModelSetup{T<:AbstractFloat,Tprog<:AbstractFloat}
    parameters::Parameter
    grid::Grid{T,Tprog}
    constants::Constants{T,Tprog}
    forcing::Forcing{T}
end

struct ModelSetup{T<:AbstractFloat}
    parameters::Parameter
    grid::Grid{T}
    constants::Constants{T}
    forcing::Forcing{T}
end

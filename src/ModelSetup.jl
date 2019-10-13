struct ModelSetup{T<:AbstractFloat}
    parameter::Parameter
    grid::Grid{T}
    constants::Constants{T}
    forcing::Forcing{T}
end

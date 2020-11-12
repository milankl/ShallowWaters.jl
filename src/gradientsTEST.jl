"""Calculates the 2nd order centred gradient in x-direction on any grid (u,v,T or q).
The size of dudx must be m-1,n compared to m,n = size(u)"""
function ∂x!(dudx::Array{T,2},u::Array{T,2}) where {T<:AbstractFloat}
    m,n = size(dudx)
    @boundscheck (m+1,n) == size(u) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            dudx[i,j] = u[i+1,j] - u[i,j]
        end
    end
end

"""Calculates the 2nd order centred gradient in y-direction on any grid (u,v,T or q).
The size of dudy must be m,n-1 compared to m,n = size(u)."""
function ∂y!(dudy::Array{T,2},u::Array{T,2}) where {T<:AbstractFloat}
    m,n = size(dudy)
    @boundscheck (m,n+1) == size(u) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            dudy[i,j] = u[i,j+1] - u[i,j]
        end
    end
end

""" ∇² is the 2nd order centred Laplace-operator ∂/∂x^2 + ∂/∂y^2.
The 1/Δ²-factor is omitted and moved into the viscosity coefficient."""
function ∇²!(du::Array{T,2},u::Array{T,2}) where {T<:AbstractFloat}
    m, n = size(du)
    @boundscheck (m+2,n+2) == size(u) || throw(BoundsError())

    minus_4 = T(-4.0)

    @inbounds for i ∈ 1:n
        for j ∈ 1:m
            du[j,i] = minus_4*u[j+1,i+1] + u[j,i+1] + u[j+2,i+1] + u[j+1,i] + u[j+1,i+2]
        end
    end
end

"""∂x is the 2nd order centred Gradient-operator ∂/∂x with grid spacing Δ (default 1)."""
function ∂x(u::Array{T,2},Δx::Real) where {T<:AbstractFloat}

    m,n = size(u)

    dudx = Array{T,2}(undef,m-1,n)
    one_over_dx = T(1.0/Δx)


    @inbounds for j ∈ 1:n
        for i ∈ 1:m-1
            dudx[i,j] = one_over_dx*(u[i+1,j] - u[i,j])
        end
    end

    return dudx
end

"""∂y is the 2nd order centred Gradient-operator ∂/∂y with grid spacing Δ (default 1)."""
function ∂y(u::Array{T,2},Δy::Real=1) where {T<:AbstractFloat}

    m,n = size(u)

    dudy = Array{T,2}(undef,m,n-1)
    one_over_dy = T(1.0/Δy)

    @inbounds for j ∈ 1:n-1
        for i ∈ 1:m
            dudy[i,j] = one_over_dy*(u[i,j+1] - u[i,j])
        end
    end

    return dudy
end

""" ∇² is the 2nd order centred Laplace-operator ∂/∂x^2 + ∂/∂y^2 with grid spacing Δ (default 1)."""
function ∇²(u::Array{T,2},Δ::Real=1) where {T<:AbstractFloat}

    m, n = size(u)
    du = Array{T,2}(undef,m-2,n-2)

    minus_4 = T(-4.0)
    one_over_dx² = T(1/Δ^2)

    @inbounds for j ∈ 2:n-1
        for i ∈ 2:m-1
            du[i-1,j-1] = one_over_dx²*(minus_4*u[i,j] + u[i,j-1] + u[i,j+1] + u[i-1,j] + u[i+1,j])
        end
    end
    return du
end

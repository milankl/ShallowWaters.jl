"""Linear interpolation of a variable u in the x-direction.
m,n = size(ux) must be m+1,n = size(u)."""
function Ix!(ux::Array{T,2},u::Array{T,2}) where {T<:AbstractFloat}
    m, n = size(ux)
    @boundscheck (m+1,n) == size(u) || throw(BoundsError())

    one_half = T(0.5)

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            ux[i,j] = one_half*(u[i+1,j] + u[i,j])
        end
    end
end

""" Linear interpolation a variable u in the y-direction.
    m,n = size(uy) must be m,n+1 = size(u)."""
function Iy!(uy::Array{T,2},u::Array{T,2}) where {T<:AbstractFloat}
    m,n = size(uy)
    @boundscheck (m,n+1) == size(u) || throw(BoundsError())

    one_half = T(0.5)

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            uy[i,j] = one_half*(u[i,j+1] + u[i,j])
        end
    end
end

""" Bilinear interpolation a variable u in x and y-direction.
m,n = size(uxy) must be m+1,n+1 = size(u). """
function Ixy!(uxy::Array{T,2},u::Array{T,2}) where {T<:AbstractFloat}
    m,n = size(uxy)
    @boundscheck (m+1,n+1) == size(u) || throw(BoundsError())

    one_quarter = T(0.25)

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            uxy[i,j] = one_quarter*(u[i,j] + u[i+1,j] + u[i,j+1] + u[i+1,j+1])
        end
    end
end

function Ix(u::Array{T,2}) where {T<:AbstractFloat}

    m,n = size(u)
    ux = Array{T,2}(undef,m-1,n)

    one_half = T(0.5)

    @inbounds for j ∈ 1:n
        for i ∈ 1:m-1
            ux[i,j] = one_half*(u[i+1,j] + u[i,j])
        end
    end

    return ux
end

function Iy(u::Array{T,2}) where {T<:AbstractFloat}

    m,n = size(u)
    uy = Array{T,2}(undef,m,n-1)

    one_half = T(0.5)

    @inbounds for j ∈ 1:n-1
        for i ∈ 1:m
            uy[i,j] = one_half*(u[i,j+1] + u[i,j])
        end
    end

    return uy
end

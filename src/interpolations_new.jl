function Ix_loop!(ux::AbstractMatrix,u::AbstractMatrix)
    #= Interpolates a variable u from any grid (u,v,T,q) to a variable ux in the x-direction.
    m,n = size(ux) must be m+1,n = size(u). =#

    m, n = size(ux)
    @boundscheck (m+1,n) == size(u) || throw(BoundsError())

    @inbounds for i ∈ 1:n
        for j ∈ 1:m
            ux[j,i] .= one_half*(u[j+1,i] .+ u[j,i])
        end
    end
end

function Ix!(ux::AbstractMatrix,u::AbstractMatrix)
    #= Interpolates a variable u from any grid (u,v,T,q) to a variable ux in the x-direction.
    m,n = size(ux) must be m+1,n = size(u). =#

    m, n = size(ux)
    @boundscheck (m+1,n) == size(u) || throw(BoundsError())

    @inbounds @views ux[:,:] .= one_half*(u[2:end,:] .+ u[1:end-1,:])
end

function Iy!(uy::AbstractMatrix,u::AbstractMatrix)
    #= Interpolates a variable u from any grid (u,v,T,q) to a variable uy in the y-direction.
    m,n = size(uy) must be m,n+1 = size(u). =#

    m, n = size(uy)
    @boundscheck (m,n+1) == size(u) || throw(BoundsError())

    @inbounds @views uy[:,:] .= one_half*(u[:,2:end] .+ u[:,1:end-1])
end

function Ixy!(uxy::AbstractMatrix,u::AbstractMatrix)
    #= Interpolates a variable u from any grid (u,v,T,q) to a variable uxy in x and y-direction.
    m,n = size(uxy) must be m+1,n+1 = size(u). =#

    m, n = size(uxy)
    @boundscheck (m+1,n+1) == size(u) || throw(BoundsError())

    @inbounds @views uxy[:,:] .= one_quarter*(u[1:end-1,1:end-1] .+ u[1:end-1,2:end] .+ u[2:end,1:end-1] .+ u[2:end,2:end])
end

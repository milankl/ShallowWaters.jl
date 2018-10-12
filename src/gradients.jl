"""Calculates the 2nd order centred gradient in x-direction on any grid (u,v,T or q).
The size of dudx must be m-1,n compared to m,n = size(u)"""
function ∂x!(dudx::AbstractMatrix,u::AbstractMatrix)
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
function ∂y!(dudy::AbstractMatrix,u::AbstractMatrix)
    m,n = size(dudy)
    @boundscheck (m,n+1) == size(u) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            dudy[i,j] = u[i,j+1] - u[i,j]
        end
    end
end

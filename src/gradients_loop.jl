function ∂x!(dudx::AbstractMatrix,u::AbstractMatrix)
    #= Calculates the gradient in x-direction on any grid (u,v,T or q).
    The size of dudx must be m-1,n compared to m,n = size(u).=#

    m, n = size(dudx)
    @boundscheck (m+1,n) == size(u) || throw(BoundsError())

    @inbounds for i ∈ 1:n
        for j ∈ 1:m
            dudx[j,i] = u[j+1,i] - u[j,i]
        end
    end
end

function ∂y!(dudy::AbstractMatrix,u::AbstractMatrix)
    #= Calculates the gradient in y-direction on any grid (u,v,T or q).
    The size of dudy must be m,n-1 compared to m,n = size(u).=#

    m, n = size(dudy)
    @boundscheck (m,n+1) == size(u) || throw(BoundsError())

    @inbounds for i ∈ 1:n
        for j ∈ 1:m
            dudx[j,i] = u[j,i+1] - u[j,i]
        end
    end
end

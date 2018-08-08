function ∇²!(du::AbstractMatrix,u::AbstractMatrix)
    #= ∇² is the Laplace-operator d/dx^2 + d/dy^2.
    The 1/dx^2 factor is omitted and moved into the viscosity coefficient. =#

    m, n = size(du)
    @boundscheck (m+2,n+2) == size(u) || throw(BoundsError())

    @views @inbounds du[:,:] .= minus_4*u[2:end-1,2:end-1] .+ u[3:end,2:end-1] .+ u[2:end-1,3:end] .+ u[1:end-2,2:end-1] .+ u[2:end-1,1:end-2]
end

# function ∇²_loop!(du::AbstractMatrix,u::AbstractMatrix)
#     #= ∇² is the Laplace-operator d/dx^2 + d/dy^2.
#     The 1/dx^2 factor is omitted and moved into the viscosity coefficient. =#
#
#     m, n = size(du)
#     @boundscheck (m+2,n+2) == size(u) || throw(BoundsError())
#
#     @inbounds for i ∈ 1:n
#         for j ∈ 1:m
#             du[j,i] = minus_4*u[j+1,i+1] + u[j,i+1] + u[j+2,i+1] + u[j+1,i] + u[j+1,i+2]
#         end
#     end
# end

#= This script defines the function in order to interpolate variables from
one grid (u,v,T,q) to another (u,v,T,q), due to the staggered Arakawa C-grid.

    The notation is Iab where a ∈ {u,v,T,q} and b ∈ {u,v,T,q}, with a being
    the departure grid and b the destination grid.

This interpolation is a 2-point averaging with coefficients (1/2,1/2) for

    ITu,ITv,IuT,IvT,Iqu,Iqv,Ivq,Iuq

and a 4-point averaging with coefficients (1/4,1/4,1/4,1/4) for

    ITq,IqT,Iuv,Ivu.

Each function takes the return variable as first argument and the variable to be
interpolated as second argument. =#

# ------------- 2-point interpolation between T and u or v. -----------------

function ITu_nonperiodic!(u::AbstractMatrix,T::AbstractMatrix)
    #= Interpolates a variable T from the T-grid to a variable u on the u-grid. =#
    @views u[:,:] .= one_half*(T[2:end,:] .+ T[1:end-1,:])
end

function ITu_periodic!(u::AbstractMatrix,T::AbstractMatrix)
    #= Interpolates a variable T from the T-grid to a variable u on the u-grid. =#
    @views u[2:end,:] .= one_half*(T[2:end,:] .+ T[1:end-1,:])
    @views u[1,:] .= one_half*(T[1,:] .+ T[end,:])
end

function ITv!(v::AbstractMatrix,T::AbstractMatrix)
    #= Interpolates a variable T from the T-grid to a variable v on the v-grid. =#
    @views v[:,:] .= one_half*(T[:,2:end] .+ T[:,1:end-1])
end

function IuT_nonperiodic!(T::AbstractMatrix,u::AbstractMatrix)
    #= Interpolates a variable u from the u-grid to a variable T on the T-grid. =#
    @views T[2:end-1,:] .= one_half * (u[2:end,:] .+ u[1:end-1,:])
    @views T[1,:] .= one_half * u[1,:]        # +0, kinematic boundary condition
    @views T[end,:] .= one_half * u[end,:]    # +0, kinematic bc
end

function IuT_periodic!(T::AbstractMatrix,u::AbstractMatrix)
    #= Interpolates a variable u from the u-grid to a variable T on the T-grid. =#
    @views T[1:end-1,:] .= one_half*(u[2:end,:] .+ u[1:end-1,:])
    @views T[end,:] .= one_half*(u[1,:]+u[end,:])
end

function IvT!(T::AbstractMatrix,v::AbstractMatrix)
    #= Interpolates a variable v from the v-grid to a variable T on the T-grid. =#
    @views T[:,2:end-1] .= one_half * (v[:,2:end] .+ v[:,1:end-1])
    @views T[:,1] .= one_half * v[:,1]        # +0, kinematic bc
    @views T[:,end] .= one_half * v[:,end]    # +0, kinematic bc
end

# ------------- 2-point interpolation between q and u or v. -----------------

function Iqu_nonperiodic!(u::AbstractMatrix,q::AbstractMatrix)
    #= Interpolates a variable q from the q-grid to a variable u on the u-grid. =#
    @views u[:,:] .= one_half * (q[2:end-1,1:end-1] .+ q[2:end-1,2:end])
end

function Iqu_periodic!(u::AbstractMatrix,q::AbstractMatrix)
    #= Interpolates a variable q from the q-grid to a variable u on the u-grid. =#
    @views u[:,:] .= one_half * (q[:,1:end-1] .+ q[:,2:end])
end

function Iqv_nonperiodic!(v::AbstractMatrix,q::AbstractMatrix)
    #= Interpolates a variable q from the q-grid to a variable v on the v-grid. =#
    @views v[:,:] .= one_half*(q[1:end-1,2:end-1] .+ q[2:end,2:end-1])
end

function Iqv_periodic!(v::AbstractMatrix,q::AbstractMatrix)
    #= Interpolates a variable q from the q-grid to a variable v on the v-grid. =#
    @views v[1:end-1,:] .= one_half * (q[1:end-1,2:end-1] .+ q[2:end,2:end-1])
    @views v[end,:] .= one_half * (q[1,2:end-1] .+ q[end,2:end-1])
end

function Ivq_nonperiodic!(q::AbstractMatrix,v::AbstractMatrix)
    #= Interpolates a variable v from the v-grid to a variable q on the q-grid. =#
    @views q[2:end-1,2:end-1] .= one_half * (v[1:end-1,:] .+ v[2:end,:])
    @views q[1,2:end-1] .= one_minus_α_half * v[1,:]
    @views q[end,2:end-1] .= one_minus_α_half * v[end,:]
    @views q[:,1] .= zeero      # if initialised with zero these lines are redundant
    @views q[:,end] .= zeero
end

function Ivq_periodic!(q::AbstractMatrix,v::AbstractMatrix)
    #= Interpolates a variable v from the v-grid to a variable q on the q-grid. =#
    @views q[2:end,2:end-1] .= one_half * (v[1:end-1,:] .+ v[2:end,:])
    @views q[1,2:end-1] .= one_half * (v[1,:] .+ v[end,:])
    @views q[:,1] .= zeero      # if initialised with zero these lines are redundant
    @views q[:,end] .= zeero
end

function Iuq_nonperiodic!(q::AbstractMatrix,u::AbstractMatrix)
    #= Interpolates a variable u from the u-grid to a variable q on the q-grid. =#
    @views q[2:end-1,2:end-1] .= one_half * (u[:,1:end-1] .+ u[:,2:end])
    @views q[1,:] .= zeero     # if initialised with zero these lines are redundant
    @views q[end,:] .= zeero
    @views q[2:end-1,1] .= one_minus_α_half * u[:,1]
    @views q[2:end-1,end] .= one_minus_α_half * u[:,end]
end

function Iuq_periodic!(q::AbstractMatrix,u::AbstractMatrix)
    #= Interpolates a variable u from the u-grid to a variable q on the q-grid. =#
    @views q[:,2:end-1] .= one_half * (u[:,1:end-1] .+ u[:,2:end])
    @views q[:,1] .= one_minus_α_half * u[:,1]
    @views q[:,end] .= one_minus_α_half * u[:,end]
end

# ------------- 4-point interpolation between (T and q) and (u and v). ---------

function ITq_nonperiodic!(q::AbstractMatrix,T::AbstractMatrix)
    #= Interpolates a variable T from the T-grid to a variable q on the q-grid. =#
    @views q[2:end-1,2:end-1] .= one_quarter * (T[1:end-1,1:end-1] .+ T[1:end-1,2:end] .+ T[2:end,1:end-1] .+ T[2:end,2:end])
    @views q[1,2:end-1] .= one_half * (T[1,1:end-1] .+ T[1,2:end])
    @views q[end,2:end-1] .= one_half * (T[end,1:end-1] .+ T[end,2:end])
    @views q[2:end-1,1] .= one_half * (T[1:end-1,1] .+ T[2:end,1])
    @views q[2:end-1,end] .= one_half * (T[1:end-1,end] .+ T[2:end,end])
    @views q[1,1] .= T[1,1]         # no gradients of T across the boundaries
    @views q[1,end] .= T[1,end]
    @views q[end,1] .= T[end,1]
    @views q[end,end] .= T[end,end]
end

function ITq_periodic!(q::AbstractMatrix,T::AbstractMatrix)
    #= Interpolates a variable T from the T-grid to a variable q on the q-grid. =#
    @views q[2:end,2:end-1] .= one_quarter * (T[1:end-1,1:end-1] .+ T[1:end-1,2:end] .+ T[2:end,1:end-1] .+ T[2:end,2:end])
    @views q[1,2:end-1] .= one_quarter * (T[1,1:end-1] .+ T[1,2:end] .+ T[end,1:end-1] .+ T[end,2:end])
    @views q[2:end,1] .= one_half * (T[1:end-1,1] .+ T[2:end,1])
    @views q[2:end,end] .= one_half * (T[1:end-1,end] .+ T[2:end,end])
    @views q[1,1] .= one_half * (T[1,1] .+ T[end,1])
    @views q[1,end] .= one_half * (T[1,end] .+ T[end,end])
end

function IqT_nonperiodic!(T::AbstractMatrix,q::AbstractMatrix)
    #= Interpolates a variable q from the q-grid to a variable T on the T-grid. =#
    @views T[:,:] .= one_quarter * (q[1:end-1,1:end-1] .+ q[1:end-1,2:end] .+ q[2:end,1:end-1] .+ q[2:end,2:end])
end

function IqT_periodic!(T::AbstractMatrix,q::AbstractMatrix)
    #= Interpolates a variable q from the q-grid to a variable T on the T-grid. =#
    @views T[1:end-1,:] .= one_quarter * (q[1:end-1,1:end-1] .+ q[1:end-1,2:end] .+ q[2:end,1:end-1] .+ q[2:end,2:end])
    @views T[end,:] .= one_quarter * (q[end,1:end-1] .+ q[end,2:end] .+ q[1,1:end-1] .+ q[1,2:end])
end

function Iuv_nonperiodic!(v::AbstractMatrix,u::AbstractMatrix)
    #= Interpolates a variable u from the u-grid to a variable v on the v-grid. =#
    @views v[2:end-1,:] .= one_quarter * (u[1:end-1,1:end-1] .+ u[1:end-1,2:end] .+ u[2:end,1:end-1] .+ u[2:end,2:end])
    @views v[1,:] .= one_quarter * (u[1,1:end-1] .+ u[1,2:end])
    @views v[end,:] .= one_quarter * (u[end,1:end-1] .+ u[end,2:end])
end

function Iuv_periodic!(v::AbstractMatrix,u::AbstractMatrix)
    #= Interpolates a variable u from the u-grid to a variable v on the v-grid. =#
    @views v[1:end-1,:] .= one_quarter * (u[1:end-1,1:end-1] .+ u[1:end-1,2:end] .+ u[2:end,1:end-1] .+ u[2:end,2:end])
    @views v[end,:] .= one_quarter * (u[1,1:end-1] .+ u[1,2:end] .+ u[end,1:end-1] .+ u[end,2:end])
end

function Ivu_nonperiodic!(u::AbstractMatrix,v::AbstractMatrix)
    #= Interpolates a variable v from the v-grid to a variable u on the u-grid. =#
    @views u[:,2:end-1] .= one_quarter * (v[1:end-1,1:end-1] .+ v[1:end-1,2:end] .+ v[2:end,1:end-1] .+ v[2:end,2:end])
    @views u[:,1] .= one_quarter * (v[1:end-1,1] .+ v[2:end,1])
    @views u[:,end] .= one_quarter * (v[1:end-1,end] .+ v[2:end,end])
end

function Ivu_periodic!(u::AbstractMatrix,v::AbstractMatrix)
    #= Interpolates a variable u from the u-grid to a variable v on the v-grid. =#
    @views u[2:end,2:end-1] .= one_quarter * (v[1:end-1,1:end-1] .+ v[1:end-1,2:end] .+ v[2:end,1:end-1] .+ v[2:end,2:end])
    @views u[1,2:end-1] .= one_quarter * (v[1,1:end-1] .+ v[1,2:end] .+ v[end,1:end-1] .+ v[end,2:end])
    @views u[2:end,1] .= one_quarter * (v[1:end-1,1] .+ v[2:end,1])
    @views u[2:end,end] .= one_quarter * (v[1:end-1,end] .+ v[2:end,end])
    @views u[1,1] .= one_quarter * (v[1,1] .+ v[end,1])
    @views u[1,end] .= one_quarter * (v[1,end] .+ v[end,end])
end

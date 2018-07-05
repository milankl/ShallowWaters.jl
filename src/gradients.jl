#= This script defines the function in order to compute the finite difference
gradient for variables on either of the grids u,v,T or q.

    The notation is Gab where a ∈ {u,v,T,q} and b ∈ {x,y}, with a being
    the departure grid and b the direction.

Each function takes the return variable as first argument and the variable to be
differentiated as second argument.

x is assumed to be the first dimension, y the second.
Both are increasing for increasing index.

α is the lateral boundary condition parameter - called lbc in its Float64 version.

    α = 0       free-slip
    0 < α < 2   partial-slip
    α = 2       no-slip

Some operators are also set-up in a non-dimensional version, i.e. dx = 1,
called _nd at the end of the function name.
=#

function GTx_nonperiodic(dTdx,T)
    #= Calculates the gradient in x-direction on the T-grid.
    The result dTdx sits on the u-grid. =#

    dTdx[:,:] = one_over_dx*(T[2:end,:] - T[1:end-1,:])
    return dTdx
end

function GTx_periodic(dTdx,T)
    #= Calculates the gradient in x-direction on the T-grid.
    The result dTdx sits on the u-grid. =#

    dTdx[2:end,:] = one_over_dx*(T[2:end,:] - T[1:end-1,:])
    dTdx[1,:] = one_over_dx*(T[1,:] - T[end,:])
    return dTdx
end

#TODO write all functions like this
#function GTx_nonperiodic_nd(dTdx::Matrix{Numtype},T::Matrix{Numtype})
#    #= Calculates the gradient in x-direction on the T-grid.
#    The result dTdx sits on the u-grid. =#
#
#    @inbounds dTdx[:,:] = T[2:end,:] - T[1:end-1,:]
# end

function GTx_nonperiodic_nd(dTdx,T)
    #= Calculates the gradient in x-direction on the T-grid.
    The result dTdx sits on the u-grid. =#

    dTdx[:,:] = T[2:end,:] - T[1:end-1,:]
    return dTdx
end

function GTx_periodic_nd(dTdx,T)
    #= Calculates the gradient in x-direction on the T-grid.
    The result dTdx sits on the u-grid. =#

    dTdx[2:end,:] = T[2:end,:] - T[1:end-1,:]
    dTdx[1,:] = T[1,:] - T[end,:]
    return dTdx
end

function GTy(dTdy,T)
    #= Calculates the gradient in y-direction on the T-grid.
    The result dTdy sits on the v-grid. =#

    dTdy[:,:] = one_over_dx*(T[:,2:end] - T[:,1:end-1])
    return dTdy
end

function GTy_nd(dTdy,T)
    #= Calculates the gradient in y-direction on the T-grid.
    The result dTdy sits on the v-grid. =#

    dTdy[:,:] = T[:,2:end] - T[:,1:end-1]
    return dTdy
end

function Gux_nonperiodic(dudx,u)
    #= Calculates the gradient in x-direction on the u-grid.
    The result dudx sits on the T-grid. =#

    dudx[2:end-1,:] = one_over_dx*(u[2:end,:] - u[1:end-1,:])
    dudx[1,:] = one_over_dx*u[1,:]              # -0, kinematic boundary condition
    dudx[end,:] = (-one_over_dx)*u[end,:]       # +0, kinematic bc
    return dudx
end

function Gux_periodic(dudx,u)
    #= Calculates the gradient in x-direction on the u-grid.
    The result dudx sits on the T-grid. =#

    dudx[1:end-1,:] = one_over_dx*(u[2:end,:]-u[1:end-1,:])
    dudx[end,:] = one_over_dx*(u[1,:]-u[end,:])
    return dudx
end

function Gux_nonperiodic_nd(dudx,u)
    #= Calculates the gradient in x-direction on the u-grid.
    The result dudx sits on the T-grid.=#

    dudx[2:end-1,:] = u[2:end,:] - u[1:end-1,:]
    dudx[1,:] = u[1,:]              # -0, kinematic boundary condition
    dudx[end,:] = u[end,:]          # +0, kinematic bc
    return dudx
end

function Gux_periodic_nd(dudx,u)
    #= Calculates the gradient in x-direction on the u-grid.
    The result dudx sits on the T-grid.=#

    dudx[1:end-1,:] = u[2:end,:]-u[1:end-1,:]
    dudx[end,:] = u[1,:]-u[end,:]
    return dudx
end

function Gvy(dvdy,v)
    #= Calculates the gradient in y-direction on the v-grid.
    The result dvdy sits on the T-grid. =#

    dvdy[:,2:end-1] = one_over_dx*(v[:,2:end] - v[:,1:end-1])
    dvdy[:,1] = one_over_dx*v[:,1]              # -0, kinematic bc
    dvdy[:,end] = (-one_over_dx)*v[:,end]       # +0, kinematic bc
    return dvdy
end

function Gvy_nd(dvdy,v)
    #= Calculates the gradient in y-direction on the v-grid.
    The result dvdy sits on the T-grid. =#

    dvdy[:,2:end-1] = v[:,2:end] - v[:,1:end-1]
    dvdy[:,1] = v[:,1]            # -0, kinematic bc
    dvdy[:,end] = -v[:,end]       # +0, kinematic bc
    return dvdy
end


function Guy_nonperiodic(dudy,u)
    #= Calculates the gradient in y-direction on the u-grid.
    The result dudy sits on the q-grid. =#

    dudy[2:end-1,2:end-1] = one_over_dx*(u[:,2:end] - u[:,1:end-1])
    dudy[2:end-1,1] = one_over_dx_α*u[:,1] #  α is the lateral boundary condition parameter
    dudy[2:end-1,end] = -one_over_dx_α*u[:,end]
    dudy[1,:] = zeero       #TODO probably redundant if initialised with zeros
    dudy[end,:] = zeero
    return dudy
end

function Guy_periodic(dudy,u)
    #= Calculates the gradient in y-direction on the u-grid.
    The result dudy sits on the q-grid. =#

    dudy[:,2:end-1] = one_over_dx*(u[:,2:end] - u[:,1:end-1])
    dudy[:,1] = one_over_dx_α*u[:,1] #  α is the lateral boundary condition parameter
    dudy[:,end] = -one_over_dx_α*u[:,end]
    return dudy
end

function Guy_nonperiodic_nd(dudy,u)
    #= Calculates the gradient in y-direction on the u-grid.
    The result dudy sits on the q-grid.
    Non-dimensional version.=#

    dudy[2:end-1,2:end-1] = u[:,2:end] - u[:,1:end-1]
    dudy[2:end-1,1] = α*u[:,1] #  α is the lateral boundary condition parameter
    dudy[2:end-1,end] = minus_α*u[:,end]
    dudy[1,:] = zeero       #TODO probably redundant if initialised with zeros
    dudy[end,:] = zeero
    return dudy
end

function Guy_periodic_nd(dudy,u)
    #= Calculates the gradient in y-direction on the u-grid.
    The result dudy sits on the q-grid.
    Non-dimensional version.=#

    dudy[:,2:end-1] = u[:,2:end] - u[:,1:end-1]
    dudy[:,1] = α*u[:,1] #  α is the lateral boundary condition parameter
    dudy[:,end] = minus_α*u[:,end]
    return dudy
end

function Gvx_nonperiodic(dvdx,v)
    #= Calculates the gradient in x-direction on the v-grid.
    The result dvdx sits on the q-grid. =#

    dvdx[2:end-1,2:end-1] = one_over_dx*(v[2:end,:] - v[1:end-1,:])
    dvdx[1,2:end-1] = one_over_dx_α*v[1,:] #  α is the lateral boundary condition parameter
    dvdx[end,2:end-1] = -one_over_dx_α*v[end,:]
    dvdx[:,1] = zeero      #TODO redundant if initialised with zeros
    dvdx[:,end] = zeero
    return dvdx
end

function Gvx_periodic(dvdx,v)
    #= Calculates the gradient in x-direction on the v-grid.
    The result dvdx sits on the q-grid. =#

    dvdx[2:end,2:end-1] = one_over_dx*(v[2:end,:] - v[1:end-1,:])
    dvdx[1,2:end-1] = one_over_dx*(v[1,:]-v[end,:])
    dvdx[:,1] = zeero      #TODO redundant if initialised with zeros
    dvdx[:,end] = zeero
    return dvdx
end

function Gvx_nonperiodic_nd(dvdx,v)
    #= Calculates the gradient in x-direction on the v-grid.
    The result dvdx sits on the q-grid. =#

    dvdx[2:end-1,2:end-1] = v[2:end,:] - v[1:end-1,:]
    dvdx[1,2:end-1] = α*v[1,:]          #  α is the lateral boundary condition parameter
    dvdx[end,2:end-1] = minus_α*v[end,:]
    dvdx[:,1] = zeero                   #TODO redundant if initialised with zeros
    dvdx[:,end] = zeero
    return dvdx
end

function Gvx_periodic_nd(dvdx,v)
    #= Calculates the gradient in x-direction on the v-grid.
    The result dvdx sits on the q-grid. =#

    dvdx[2:end,2:end-1] = v[2:end,:] - v[1:end-1,:]
    dvdx[1,2:end-1] = v[1,:]-v[end,:]
    dvdx[:,1] = zeero      #TODO redundant if initialised with zeros
    dvdx[:,end] = zeero
    return dvdx
end

function Gqx_nonperiodic(dqdx,q)
    #= Calculates the gradient in x-direction on the q-grid.
    The result dqdx sits on the v-grid. =#

    dqdx[:,:] = one_over_dx*(q[2:end,2:end-1] - q[1:end-1,2:end-1])
    return dqdx
end

function Gqx_periodic(dqdx,q)
    #= Calculates the gradient in x-direction on the q-grid.
    The result dqdx sits on the v-grid. =#

    dqdx[1:end-1,:] = one_over_dx*(q[2:end,2:end-1] - q[1:end-1,2:end-1])
    dqdx[end,:] = one_over_dx*(q[1,2:end-1] - q[end,2:end-1])
    return dqdx
end

function Gqx_nonperiodic_nd(dqdx,q)
    #= Calculates the gradient in x-direction on the q-grid.
    The result dqdx sits on the v-grid. =#

    dqdx[:,:] = q[2:end,2:end-1] - q[1:end-1,2:end-1]
    return dqdx
end

function Gqx_periodic_nd(dqdx,q)
    #= Calculates the gradient in x-direction on the q-grid.
    The result dqdx sits on the v-grid. =#

    dqdx[1:end-1,:] = q[2:end,2:end-1] - q[1:end-1,2:end-1]
    dqdx[end,:] = q[1,2:end-1] - q[end,2:end-1]
    return dqdx
end

function Gqy_nonperiodic(dqdy,q)
    #= Calculates the gradient in y-direction on the q-grid.
    The result dqdy sits on the u-grid. =#

    dqdy[:,:] = one_over_dx*(q[2:end-1,2:end] - q[2:end-1,1:end-1])
    return dqdy
end

function Gqy_periodic(dqdy,q)
    #= Calculates the gradient in y-direction on the q-grid.
    The result dqdy sits on the u-grid. =#

    dqdy[:,:] = one_over_dx*(q[:,2:end] - q[:,1:end-1])
    return dqdy
end

function Gqy_nonperiodic_nd(dqdy,q)
    #= Calculates the gradient in y-direction on the q-grid.
    The result dqdy sits on the u-grid. =#

    dqdy[:,:] = q[2:end-1,2:end] - q[2:end-1,1:end-1]
    return dqdy
end

function Gqy_periodic_nd(dqdy,q)
    #= Calculates the gradient in y-direction on the q-grid.
    The result dqdy sits on the u-grid. =#

    dqdy[:,:] = q[:,2:end] - q[:,1:end-1]
    return dqdy
end

#= This script defines the function in order to compute the finite difference
gradient for variables on either of the grids u,v,T or q.

    The notation is Gab where a ∈ {u,v,T,q} and b ∈ {x,y}, with a being
    the departure grid and b the direction.

Each function takes the return variable as first argument and the variable to be
differentiated as second argument. =#

function GTx_nonperiodic(dTdx,T)
    #= Calculates the gradient in x-direction on the T-grid.
    The result dTdx sits on the u-grid. =#

    # non-periodic version, assumes x is the first dimension
    dTdx[:,:] = one_over_dx*(T[2:end,:] - T[1:end-1,:])
    return dTdx
end

function GTx_periodic(dTdx,T)
    #= Calculates the gradient in x-direction on the T-grid.
    The result dTdx sits on the u-grid. =#

    # periodic version, assumes x is the first dimesion
    dTdx[2:end,:] = one_over_dx*(T[2:end,:] - T[1:end-1,:])
    dTdx[1,:] = one_over_dx*(T[1,:] - T[end,:])
    return dTdx
end

function GTy(dTdy,T)
    #= Calculates the gradient in y-direction on the T-grid.
    The result dTdy sits on the v-grid. =#

    # non-periodic version, y is second dimension increasing for increasing
    # index
    dTdy[:,:] = one_over_dx*(T[:,2:end] - T[:,1:end-1])
    return dTdy
end

function Gux_nonperiodic(dudx,u)
    #= Calculates the gradient in x-direction on the u-grid.
    The result dudx sits on the T-grid. =#

    dudx[2:end-1,:] = one_over_dx*(u[2:end,:]-u[1:end-1,:])
    dudx[1,:] = one_over_dx*u[1,:]    # -0, kinematic boundary condition
    dudx[end,:] = (-one_over_dx)*u[end,:] # +0, kinematic bc
    return dudx
end

function Gux_periodic(dudx,u)
    #= Calculates the gradient in x-direction on the u-grid.
    The result dudx sits on the T-grid. =#

    dudx[1:end-1,:] = one_over_dx*(u[2:end,:]-u[1:end-1,:])
    dudx[end,:] = one_over_dx*(u[end,:]-u[1,:])
    return dudx
end

function Gvy(dvdy,v)
    #= Calculates the gradient in y-direction on the v-grid.
    The result dvdy sits on the T-grid. =#

    dvdy[:,2:end-1] = one_over_dx*(v[:,2:end] - v[:,1:end-1])
    dvdy[:,1] = one_over_dx*v[:,1] # -0, kinematic bc
    dvdy[:,end] = (-one_over_dx)*v[:,end] # +0, kinematic bc
    return dvdy
end

function Guy(dudy,u)
    #= Calculates the gradient in y-direction on the u-grid.
    The result dudy sits on the q-grid. =#

    dudy[:,2:end-1] = one_over_dx*(u[:,2:end] - u[:,1:end-1])
    dudy[:,1] = (one_over_dx*α)*u[:,1] #  α is the lateral boundary condition parameter
    dudy[:,end] = (-one_over_dx*α)*u[:,end]
    return dudy
end

function Gvx_nonperiodic(dvdx,v)
    #= Calculates the gradient in x-direction on the v-grid.
    The result dvdx sits on the q-grid. =#

    dvdx[2:end-1,:] = one_over_dx*(v[2:end,:] - v[1:end-1,:])
    dvdx[1,:] = (one_over_dx*α)*v[1,:] #  α is the lateral boundary condition parameter
    dvdx[end,:] = (-one_over_dx*α)*v[end,:]
    return dvdx
end



function Gvx_periodic(dvdx,v)
    #= Calculates the gradient in x-direction on the v-grid.
    The result dvdx sits on the q-grid. =#

    dvdx[2:end,:] = one_over_dx*(v[2:end,:] - v[1:end-1,:])
    dvdx[1,:] = one_over_dx*(v[1,:]-v[end,:])
    return dvdx
end

function Gqx_nonperiodic(dqdx,q)
    #= Calculates the gradient in x-direction on the q-grid.
    The result dqdx sits on the v-grid. =#

    dqdx[:,:] = one_over_dx*(q[2:end,:] - q[1:end-1,:])
    return dqdx
end

function Gqx_periodic(dqdx,q)
    #= Calculates the gradient in x-direction on the q-grid.
    The result dqdx sits on the v-grid. =#

    dqdx[1:end-1,:] = one_over_dx*(q[2:end,:] - q[1:end-1,:])
    dqdx[end,:] = one_over_dx*(q[1,:] - q[end,:])
    return dqdx
end

function Gqy(dqdy,q)
    #= Calculates the gradient in y-direction on the q-grid.
    The result dqdy sits on the u-grid. =#

    dqdy[:,:] = one_over_dx*(q[:,2:end] - q[:,1:end-1])    
    return dqdy
end

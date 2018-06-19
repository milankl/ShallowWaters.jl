#= This script defines the function in order to compute the finite difference
gradient for variables on either of the grids u,v,T or q.

    The notation is Gab where a ∈ {u,v,T,q} and b ∈ {x,y}, with a being
    the departure grid and b the direction.

Each function takes the return variable as first argument and the variable to be
differentiated as second argument. =#

#TODO BOUNDARY CONDITIONS

function GTx(dTdx,T)
    #= Calculates the gradient in x-direction on the T-grid.
    The result dTdx sits on the u-grid. =#
    return dTdx
end

function GTy(dTdy,T)
    #= Calculates the gradient in y-direction on the T-grid.
    The result dTdy sits on the v-grid. =#
    return dTdy
end

function Gux(dudx,u)
    #= Calculates the gradient in x-direction on the u-grid.
    The result dudx sits on the T-grid. =#
    return dudx
end

function Guy(dudy,u)
    #= Calculates the gradient in y-direction on the u-grid.
    The result dudy sits on the q-grid. =#
    return dudy
end

function Gvx(dvdx,v)
    #= Calculates the gradient in x-direction on the v-grid.
    The result dvdx sits on the q-grid. =#
    return dvdx
end

function Gvy(dvdy,v)
    #= Calculates the gradient in y-direction on the v-grid.
    The result dvdy sits on the T-grid. =#
    return dvdy
end

function Gqx(dqdx,q)
    #= Calculates the gradient in x-direction on the q-grid.
    The result dqdx sits on the v-grid. =#
    return dqdx
end

function Gqy(dqdy,q)
    #= Calculates the gradient in y-direction on the q-grid.
    The result dqdy sits on the v-grid. =#
    return dqdy
end

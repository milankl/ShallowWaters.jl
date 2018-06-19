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

function ITu(u,T)
    #= Interpolates a variable T from the T-grid to a variable u on the u-grid. =#
    return u
end

function ITv(v,T)
    #= Interpolates a variable T from the T-grid to a variable v on the v-grid. =#
    return v
end

function IuT(T,u)
    #= Interpolates a variable u from the u-grid to a variable T on the T-grid. =#
    return T
end

function IvT(T,v)
    #= Interpolates a variable v from the v-grid to a variable T on the T-grid. =#
    return T
end

# ------------- 2-point interpolation between q and u or v. -----------------

function Iqu(u,q)
    #= Interpolates a variable q from the q-grid to a variable u on the u-grid. =#
    return u
end

function Iqv(v,q)
    #= Interpolates a variable q from the q-grid to a variable v on the v-grid. =#
    return v
end

function Ivq(q,v)
    #= Interpolates a variable v from the v-grid to a variable q on the q-grid. =#
    return q
end

function Iuq(q,u)
    #= Interpolates a variable u from the u-grid to a variable q on the q-grid. =#
    return q
end

# ------------- 4-point interpolation between (T and q) and (u and v). ---------

function ITq(q,T)
    #= Interpolates a variable T from the T-grid to a variable q on the q-grid. =#
    return q
end

function IqT(T,q)
    #= Interpolates a variable q from the q-grid to a variable T on the T-grid. =#
    return T

function Iuv(v,u)
    #= Interpolates a variable u from the u-grid to a variable v on the v-grid. =#
    return v
end

function Ivu(u,v)
    #= Interpolates a variable v from the v-grid to a variable u on the u-grid. =#
    return u
end

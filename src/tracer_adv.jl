"""Computes the departure point for semi-Lagrangian advection following Diamantakis, 2014.
u,v are assumed to be the time averaged velocities over the previous advection time step.
(Presumably need to be changed to 2nd order extrapolation in case the tracer is not passive)

Uses fixed-point iteration (sofar only once) to find the departure point."""
function departure!(u,v,u_T,v_T,um,vm,um_T,vm_T,uinterp,vinterp,xd,yd)
    # u,v is t + dtadv, um,vm are averaged over (t,t_dtadv)

    # interpolate u,um,v,vm onto the T-grid
    Ix!(u_T,u)
    Ix!(um_T,um)
    Iy!(v_T,v)
    Iy!(vm_T,v)

    # initial guess - mid point
    backtraj!(xd,xxT,one_half*dtadv,u_T)
    backtraj!(yd,yyT,one_half*dtadv,v_T)

    # interpolate um,vm onto mid-point
    interp!(uinterp,um_T,xd,yd)
    interp!(vinterp,vm_T,xd,yd)

    # update departure point
    backtraj!(xd,xxT,dtadv,uinterp)
    backtraj!(yd,yyT,dtadv,vinterp)
end

""" Solves the trajectory equation for a given arrival point ra (this can be either x or y),
a time step dt and the velocity uv (this can be u or v). All matrices have to be of the same size."""
function backtraj!(rd::AbstractMatrix,ra::AbstractMatrix,dt::Real,uv::AbstractMatrix)
    m,n = size(rd)
    @boundscheck (m,n) == size(ra) || throw(BoundsError())
    @boundscheck (m,n) == size(uv) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            rd[i,j] = ra[i,j] - dt*uv[i,j]
        end
    end
end

function interp!(uinterp,um,xdm,ydm)

end

function adv_sst!(ssti,sst,xd,yd)

end

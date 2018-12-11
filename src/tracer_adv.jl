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
    Iy!(vm_T,vm)

    # initial guess - mid point
    backtraj_ump!(xd,xxT,one_half*dtadvu,u_T)
    backtraj_vmp!(yd,yyT,one_half*dtadvv,v_T)
    #backtraj_ump!(xd,xxT,dtadvu,um_T)
    #backtraj_vmp!(yd,yyT,dtadvv,vm_T)

    # # interpolate um,vm onto mid-point
    # #TODO make one function currently not possible because of different matrix sizes
    interp_u!(uinterp,um_T,xd,yd)
    interp_v!(vinterp,vm_T,xd,yd)
    #
    # # update departure point
    backtraj!(xd,xxT,dtadvu,uinterp)
    backtraj!(yd,yyT,dtadvv,vinterp)
end

""" Solves the trajectory equation for a given arrival point xa, a time step dt and the velocity u.
xd, xa sit on the T-grid without halo, u is interpolated from u-grid (with halo) onto T-grd.
Respective halo points will be ignored via indexing."""
function backtraj_ump!(xd::AbstractMatrix,xa::AbstractMatrix,dt::Real,u::AbstractMatrix)
    m,n = size(xd)
    @boundscheck (m,n) == size(xa) || throw(BoundsError())
    @boundscheck (m+2+ep,n+4) == size(u) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            xd[i,j] = xa[i,j] - dt*u[i+1+ep,j+2]
        end
    end
end

""" Solves the trajectory equation for a given arrival point ya, a time step dt and the velocity v.
yd, ya sit on the T-grid without halo, v is interpolated from v-grid (with halo) onto T-grid.
Respective halo points will be ignored via indexing."""
function backtraj_vmp!(yd::AbstractMatrix,ya::AbstractMatrix,dt::Real,v::AbstractMatrix)
    m,n = size(yd)
    @boundscheck (m,n) == size(ya) || throw(BoundsError())
    @boundscheck (m+4,n+2) == size(v) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            yd[i,j] = ya[i,j] - dt*v[i+2,j+1]
        end
    end
end

""" Solves the trajectory equation for a given arrival point ra (this can be either x or y),
a time step dt and the velocity uv (this can be u or v). All matrices have to be of the same size.
Update version, where all matrices have same size."""
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

""" Interpolates the matrix u into the matrix ui, where xx, yy specify the coordinates as indices (including fraction.
Interpolation only onto the inner entries of ui. (They will be copied back later via the ghostpoint function)"""
function interp_u!(ui,u,xx,yy)
    m,n = size(ui)
    @boundscheck (m+2+ep,n+4) == size(u) || throw(BoundsError())
    @boundscheck (m,n) == size(xx) || throw(BoundsError())
    @boundscheck (m,n) == size(yy) || throw(BoundsError())

    # clip to avoid indices beyond [1,m]x[1,n]
    clip_wrap!(xx,Numtype(-ep),Numtype(nx+1))
    clip!(yy,Numtype(-1),Numtype(ny+2))

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            # floor is not defined for posits...
            xi = Int(floor(Float64(xx[i,j])))
            yi = Int(floor(Float64(yy[i,j])))
            k = xi+1+ep   #+1 to account for the size difference of u,ui
            l = yi+2
            x0 = xx[i,j] - xi
            y0 = yy[i,j] - yi
            ui[i,j] = bilin(u[k,l],u[k+1,l],u[k,l+1],u[k+1,l+1],x0,y0)
        end
    end
end

""" Interpolates the matrix u into the matrix ui, where xx, yy specify the coordinates as indices (including fraction.
Interpolation only onto the inner entries of ui. (They will be copied back later via the ghostpoint function)"""
function interp_v!(vi,v,xx,yy)
    m,n = size(vi)
    @boundscheck (m+4,n+2) == size(v) || throw(BoundsError())
    @boundscheck (m,n) == size(xx) || throw(BoundsError())
    @boundscheck (m,n) == size(yy) || throw(BoundsError())

    # clip to avoid indices beyond [1,m]x[1,n]
    clip_wrap!(xx,Numtype(-1),Numtype(nx+2))
    clip!(yy,Numtype(0),Numtype(ny+1))

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            # floor is not defined for posits...
            xi = Int(floor(Float64(xx[i,j])))
            yi = Int(floor(Float64(yy[i,j])))
            k = xi+2   #+1 to account for the size difference of u,ui
            l = yi+1
            x0 = xx[i,j] - xi
            y0 = yy[i,j] - yi
            vi[i,j] = bilin(v[k,l],v[k+1,l],v[k,l+1],v[k+1,l+1],x0,y0)
        end
    end
end

""" At the moment this is except for matrix sizes the same as interp!, will be changed in the future?
    #TODO think about using a higher order interpolation?"""
function adv_sst!(ssti,sst,xx,yy)
    m,n = size(ssti)
    @boundscheck (m,n) == size(sst) || throw(BoundsError())
    @boundscheck (m-2*halosstx,n-2*halossty) == size(xx) || throw(BoundsError())
    @boundscheck (m-2*halosstx,n-2*halossty) == size(yy) || throw(BoundsError())

    clip_wrap!(xx,Numtype(1-halosstx),Numtype(nx+halosstx))
    clip!(yy,Numtype(1-halossty),Numtype(ny+halossty))

    for j ∈ halossty+1:n-halossty
        for i ∈ halosstx+1:m-halosstx
            xi = Int(floor(Float64(xx[i-halosstx,j-halossty])))   # departure point
            yi = Int(floor(Float64(yy[i-halosstx,j-halossty])))
            k = xi+halosstx
            l = yi+halossty
            x0 = xx[i-halosstx,j-halossty] - xi
            y0 = yy[i-halosstx,j-halossty] - yi
            ssti[i,j] = bilin(sst[k,l],sst[k+1,l],sst[k,l+1],sst[k+1,l+1],x0,y0)
        end
    end
end

"""Bilinear interpolation on (x,y) in the unit square [0,1]x[0,1].
The values at the corners are f00 = f(0,0), f01 = f(0,1), etc."""
function bilin(f00::Real,f10::Real,f01::Real,f11::Real,x::Real,y::Real)
    return f00*(oone-x)*(oone-y) + f10*x*(oone-y) + f01*(oone-x)*y + f11*x*y
end


"""Tracer relaxation."""
function tracer_relax!(sst,sst_ref)
    m,n = size(sst)
    @boundscheck (m-2*halosstx,n-2*halossty) == size(sst_ref) || throw(BoundsError())

    @inbounds for j ∈ halossty:n-halossty
        for i ∈ halosstx:m-halosstx
            sst[i,j] += r_SST*(sst_ref[i-halosstx,j-halossty] - sst[i,j])
        end
    end
end

"""Clips all values of Matrix X in the range [a,b)."""
function clip!(X::AbstractMatrix,a::Real,b::Real)
    if minimum(X) < a || maximum(X) >= b
        #println("Limits exceed matrix dimensions. Clipping...")
        X[X .< a] .= a
        X[X .>= b] .= b
    end
    return nothing
end

"""Clips all values of Matrix X in the range [a,b) with wrap-around behaviour:
x* = x + (b-a) for x < a,   and
x* = x + (a-b) for x >= b"""
function clip_wrap!(X::AbstractMatrix,a::Real,b::Real)
    if minimum(X) < a || maximum(X) >= b
        #println("Limits exceed matrix dimensions. Wrapping...")
        X[X .< a] .+= (b-a)
        X[X .>= b] .+= (a-b)
    end
    return nothing
end

#TODO for no domain decomposition the clipping function could involve a wrap-around
#TODO for periodic boundary conditions (probably unnecessary )

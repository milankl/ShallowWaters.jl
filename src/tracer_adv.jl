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
    #backtraj_ump!(xd,xxT,one_half*dtadvu,um_T)
    #backtraj_vmp!(yd,yyT,one_half*dtadvv,vm_T)
    backtraj_ump!(xd,xxT,dtadvu,um_T)
    backtraj_vmp!(yd,yyT,dtadvv,vm_T)

    # # interpolate um,vm onto mid-point
    # #TODO make one function currently not possible because of different matrix sizes
    # interp_u!(uinterp,um_T,xd,yd)
    # interp_v!(vinterp,vm_T,xd,yd)
    #
    # # update departure point
    # backtraj!(xd,xxT,dtadvu,uinterp)
    # backtraj!(yd,yyT,dtadvv,vinterp)
end

""" Solves the trajectory equation for a given arrival point xa, a time step dt and the velocity u.
xd, xa sit on the T-grid without halo, u is interpolated from u-grid (with halo) onto T-grd.
Respective halo points will be ignored via indexing."""
function backtraj_ump!(xd::AbstractMatrix,xa::AbstractMatrix,dt::Real,u::AbstractMatrix)
    m,n = size(xd)
    @boundscheck (m,n) == size(xa) || throw(BoundsError())
    @boundscheck (m+2+ep,n+4) == size(u) || throw(BoundsError())

    #@inbounds for j ∈ 1:n
    for j ∈ 1:n
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

    #@inbounds for j ∈ 1:n
    for j ∈ 1:n
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

    #@inbounds for j ∈ 1:n
    for j ∈ 1:n
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

    # max/min to avoid indices beyond [1,m]x[1,n]
    if minimum(xx) < -ep || maximum(xx) >= (nx+1)
        throw(error("Interpolation error: Index exceeds matrix dimensions."))
        # TODO think what to do then
    end

    if minimum(yy) < -1 || maximum(yy) >= (ny+2)
        throw(error("Interpolation error: Index exceeds matrix dimensions."))
        # TODO think what to do then
    end

    #@inbounds for j ∈ 1:n
    for j ∈ 1:n
        for i ∈ 1:m
            k = Int(floor(xx[i,j]))+1+ep   #+1 to account for the size difference of u,ui
            l = Int(floor(yy[i,j]))+2
            x0 = xx[i,j] % 1
            y0 = yy[i,j] % 1
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

    # max/min to avoid indices beyond [1,m]x[1,n]
    if minimum(xx) < -1 || maximum(xx) >= (nx+2)
        throw(error("Interpolation error: Index exceeds matrix dimensions."))
        # TODO think what to do then
    end

    if minimum(yy) < 0 || maximum(yy) >= (ny+1)
        throw(error("Interpolation error: Index exceeds matrix dimensions."))
        # TODO think what to do then
    end

    #@inbounds for j ∈ 1:n
    for j ∈ 1:n
        for i ∈ 1:m
            k = Int(floor(xx[i,j]))+2   #+1 to account for the size difference of u,ui
            l = Int(floor(yy[i,j]))+1
            x0 = xx[i,j] % 1
            y0 = yy[i,j] % 1
            vi[i,j] = bilin(v[k,l],v[k+1,l],v[k,l+1],v[k+1,l+1],x0,y0)
        end
    end
end

""" At the moment this is except for matrix sizes the same as interp!, will be changed in the future?
    #TODO think about using a higher order interpolation?"""
function adv_sst!(ssti,sst,xx,yy)
    m,n = size(ssti)
    @boundscheck (m,n) == size(sst) || throw(BoundsError())
    @boundscheck (m-2,n-2) == size(xx) || throw(BoundsError())
    @boundscheck (m-2,n-2) == size(yy) || throw(BoundsError())

    # max/min to avoid indices beyond [1,m]x[1,n]
    if minimum(xx) < 0 || maximum(xx) >= (nx+1)
        throw(error("Interpolation error: Index exceeds matrix dimensions."))
        # TODO think what to do then
    end

    if minimum(yy) < 0 || maximum(yy) >= (ny+1)
        throw(error("Interpolation error: Index exceeds matrix dimensions."))
        # TODO think what to do then
    end

    #@inbounds for j ∈ 2:n-1
    for j ∈ 2:n-1
        for i ∈ 2:m-1
            k = Int(floor(xx[i-1,j-1]))+1
            l = Int(floor(yy[i-1,j-1]))+1
            x0 = xx[i-1,j-1] % 1
            y0 = yy[i-1,j-1] % 1
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
    @boundscheck (m-2,n-2) == size(sst_ref) || throw(BoundsError())

    @inbounds for i ∈ 2:m-1
        for j ∈ 2:n-1
            sst[i,j] += r_SST*(sst_ref[i-1,j-1] - sst[i,j])
        end
    end
end

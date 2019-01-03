"""Computes the departure point for semi-Lagrangian advection following Diamantakis, 2014.
u,v are assumed to be the time averaged velocities over the previous advection time step.
(Presumably need to be changed to 2nd order extrapolation in case the tracer is not passive)

Uses fixed-point iteration once to find the departure point."""
function departure!(u,v,u_T,v_T,um,vm,um_T,vm_T,uinterp,vinterp,xd,yd)
    # u,v is t + dtadv, um,vm are averaged over (t,t_dtadv)

    # interpolate u,um,v,vm onto the T-grid
    Ix!(u_T,u)
    Ix!(um_T,um)
    Iy!(v_T,v)
    Iy!(vm_T,vm)

    # initial guess for departure point - mid point
    backtraj!(xd,xxT,one_half*dtadvu,u_T)
    backtraj!(yd,yyT,one_half*dtadvv,v_T)

    # interpolate um,vm onto mid-point
    interp_uv!(uinterp,um_T,xd,yd)
    interp_uv!(vinterp,vm_T,xd,yd)

    # update departure point
    backtraj!(xd,xxT,dtadvu,uinterp)
    backtraj!(yd,yyT,dtadvv,vinterp)
end

""" Solves the trajectory equation for a given arrival point ra (this can be either x or y),
a time step dt and the velocity uv (this can be u or v). One function for three cases

(i) u is interpolated from u-grid with halo onto T-grid
(ii) v is interpolated from v-grid with halo onto T-grid
(iii) u or v already on the T-grid: All matrices have same size.

Uses relative grid nodes in the departure points rd, such that actually solved is rd = 0 - dt*uv.
The indices i,j of rd determine the arrival grid point."""
function backtraj!(rd::AbstractMatrix,ra::AbstractMatrix,dt::Real,uv::AbstractMatrix)
    m,n = size(rd)
    @boundscheck (m,n) == size(ra) || throw(BoundsError())

    if (m,n) == size(uv)            # update departue point case, matrices have same size
        ishift = 0
        jshift = 0
    elseif (m+4,n+2) == size(uv)    # v-vel mid-point case, v has halo
        ishift = 2
        jshift = 1
    elseif (m+2+ep,n+4) == size(uv) # u-vel mid-point case, u has halo
        ishift = 1+ep
        jshift = 2
    else
        throw(BoundsError())
    end

    # relative grid means rd = 0 - dt*uv. The arrival information is stored in
    # the indices i,j of rd: rd[2,3] => -0.5,-0.5 means for the arrival grid node (2,3)
    # the departure is (2-0.5,3-0.5) = (1.5,2.5)
    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            rd[i,j] = -dt*uv[i+ishift,j+jshift]
        end
    end
end

""" Interpolates the matrix uv into the matrix uvi, where xx, yy specify the coordinates as indices (including fraction.
Interpolation only onto the inner entries of uvi. (They will be copied back later via the ghostpoint function).
Two cases
    (i) u velocities: from u-grid with halo to T-grid
    (ii) v velocities: from v-grid with halo to T-grid."""
function interp_uv!(uvi::AbstractMatrix,uv::AbstractMatrix,xx::AbstractMatrix,yy::AbstractMatrix)
    m,n = size(uvi)
    @boundscheck (m,n) == size(xx) || throw(BoundsError())
    @boundscheck (m,n) == size(yy) || throw(BoundsError())

    if (m+2+ep,n+4) == size(uv)      # u case
        ishift = 1+ep
        jshift = 2
        clip_x!(xx,Numtype(-ep),Numtype(nx+1))
    elseif (m+4,n+2) == size(uv)    # v case
        ishift = 2
        jshift = 1
        clip_x!(xx,Numtype(-1),Numtype(nx+2))
    else
        throw(BoundsError())
    end

    #TODO take clip_x somehow out of the if-clause?
    #clip_x!(xx,Numtype(1-ishift),Numtype(nx+2))
    clip_y!(yy,Numtype(1-jshift),Numtype(ny+jshift))

    for j ∈ 1:n
        for i ∈ 1:m
            # floor is not defined for posits...
            xi = Int(floor(Float64(xx[i,j])))   # departure point indices lower left corner within grid cell
            yi = Int(floor(Float64(yy[i,j])))
            k = xi+ishift+i       # indices of uv
            l = yi+jshift+j
            x0 = xx[i,j] - xi   # coordinates within grid cell
            y0 = yy[i,j] - yi
            uvi[i,j] = bilin(uv[k,l],uv[k+1,l],uv[k,l+1],uv[k+1,l+1],x0,y0)
        end
    end
end

""" At the moment this is except for matrix sizes the same as interp!, will be changed in the future?
    #TODO think about using a higher order interpolation?"""
function adv_sst!(ssti::AbstractMatrix,sst::AbstractMatrix,xx::AbstractMatrix,yy::AbstractMatrix)
    m,n = size(ssti)
    @boundscheck (m,n) == size(sst) || throw(BoundsError())
    @boundscheck (m-2*halosstx,n-2*halossty) == size(xx) || throw(BoundsError())
    @boundscheck (m-2*halosstx,n-2*halossty) == size(yy) || throw(BoundsError())

    clip_x!(xx,Numtype(1-halosstx),Numtype(nx+halosstx-1))
    clip_y!(yy,Numtype(1-halossty+1),Numtype(ny+halossty-1))

    for j ∈ 1:n-2*halossty
        for i ∈ 1:m-2*halosstx
            xi = Int(floor(Float64(xx[i,j])))   # departure point
            yi = Int(floor(Float64(yy[i,j])))
            k = xi+i+halosstx
            l = yi+j+halossty
            x0 = xx[i,j] - xi                   # subtraction Numtype - Int
            y0 = yy[i,j] - yi
            if x0 < 0 || x0 > 1 || y0 < 0 || y0 > 1
                println("x0,y0 out of bounds")
            end
            try
                sst[k,l]
            catch
                println((k,l,x0,y0))
            end
            ssti[i,j] = bilin(sst[k,l],sst[k+1,l],sst[k,l+1],sst[k+1,l+1],x0,y0)
        end
    end
end

"""Bilinear interpolation on (x,y) in the unit square [0,1]x[0,1].
The values at the corners are f00 = f(0,0), f01 = f(0,1), etc."""
function bilin(f00::Real,f10::Real,f01::Real,f11::Real,x::Real,y::Real)
    # oone is Numtype(1)
    return f00*(oone-x)*(oone-y) + f10*x*(oone-y) + f01*(oone-x)*y + f11*x*y
end


"""Tracer relaxation."""
function tracer_relax!(sst::AbstractMatrix,sst_ref::AbstractMatrix)
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

"""Clips relative departure points x-coordinates."""
function clip_rel_x!(xx::AbstractMatrix,a::Real,b::Real)
    m,n = size(xx)

    for j ∈ 1:n
        for i ∈ 1:m
            # this is a comparison of Numtype+Int < Numtype
            if xx[i,j]+i < a
                xx[i,j] = Numtype(a-i)
            elseif xx[i,j]+i >= b
                xx[i,j] = Numtype(b-i)
            end
        end
    end
end

"""Clips relative departure points y-coordinates."""
function clip_rel_y!(yy::AbstractMatrix,a::Real,b::Real)
    m,n = size(yy)

    for j ∈ 1:n
        for i ∈ 1:m
            # this is a comparison of Numtype+Int < Numtype
            if yy[i,j]+j < a
                yy[i,j] = Numtype(a-j)
            elseif yy[i,j]+j >= b
                yy[i,j] = Numtype(b-j)
            end
        end
    end
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

if bc_x == "periodic"
    #TODO
    clip_x! = clip_wrap!
    clip_y! = clip!
else
    clip_x! = clip_rel_x!
    clip_y! = clip_rel_y!
end

function tracer!(   i::Integer,
                    Prog::PrognosticVar,
                    Diag::DiagnosticVar,
                    S::ModelSetup)

    @unpack tracer_advection = S.parameters
    @unpack nadvstep_half,nadvstep = S.grid

    # mid point (in time) velocity for the advective time step
    if tracer_advection && ((i+nadvstep_half) % nadvstep) == 0

        @unpack u,v = Prog
        @unpack um,vm = Diag.SemiLagrange

        copyto!(um,u)
        copyto!(vm,v)
    end

    if tracer_advection && (i % nadvstep) == 0
        #departure!(u,v,u_T,v_T,um,vm,um_T,vm_T,uinterp,vinterp,xd,yd)
        departure!(Prog,Diag)
        #adv_sst!(ssti,sst,xd,yd)
        adv_sst!(Prog,Diag)

        # if tracer_relaxation
        #     tracer_relax!(ssti,sst_ref,SSTγ)
        # end
        # if tracer_consumption
        #     tracer_consumption!(ssti)
        # end

        ghost_points_sst!(ssti,S)
        copyto!(sst,ssti)

        #TODO tracer conservation?
        #println(mean(sst[halosstx+1:end-halosstx,halossty+1:end-halossty].*h[haloη+1:end-haloη,haloη+1:end-haloη]))
    end
end


"""Computes the departure point for semi-Lagrangian advection following Diamantakis, 2014.
u,v are assumed to be the time averaged velocities over the previous advection time step.
(Presumably need to be changed to 2nd order extrapolation in case the tracer is not passive)

Uses fixed-point iteration once to find the departure point."""
function departure!(Prog::PrognosticVar,
                    Diag::DiagnosticVar,
                    S::ModelSetup)

    @unpack u,v = Prog
    @unpack u_T,v_T,um,vm,um_T,vm_T,uinterp,vinterp,xd,yd = Diag.SemiLagrange
    @unpack dtadvu,dtadv,half_dtadvu,half_dtadvv = S.grid

    # u,v is t + dtadv, um,vm are averaged over (t,t_dtadv)
    # interpolate u,um,v,vm onto the T-grid
    Ix!(u_T,u)
    Ix!(um_T,um)
    Iy!(v_T,v)
    Iy!(vm_T,vm)

    # simple mid-point in time advection
    #backtraj!(xd,dtadvu,um_T)
    #backtraj!(yd,dtadvv,vm_T)

    # fixed-point iteration for departure point
    # initial guess for departure point - mid point
    backtraj!(xd,half_dtadvu,u_T,ep)
    backtraj!(yd,half_dtadvv,v_T,ep)

    # interpolate um,vm onto mid-point
    interp_uv!(uinterp,um_T,xd,yd)
    interp_uv!(vinterp,vm_T,xd,yd)

    # update departure point
    backtraj!(xd,dtadvu,uinterp,ep)
    backtraj!(yd,dtadvv,vinterp,ep)
end

""" Solves the trajectory equation for a given time step dt and the velocity uv (this can be u or v).
One function for three cases

(i) u is interpolated from u-grid with halo onto T-grid
(ii) v is interpolated from v-grid with halo onto T-grid
(iii) u or v already on the T-grid: All matrices have same size.

Uses relative grid nodes in the departure points rd, such that actually solved is rd = 0 - dt*uv.
The indices i,j of rd determine the arrival grid point."""
function backtraj!( rd::Array{T,2},
                    dt::T,
                    uv::Array{T,2},
                    ep::Int) where {T<:AbstractFloat}
    m,n = size(rd)

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
function interp_uv!(uvi::AbstractMatrix,
                    uv::AbstractMatrix,
                    xx::AbstractMatrix,
                    yy::AbstractMatrix,
                    ep::Int)

    m,n = size(uvi)
    muv,nuv = size(uv)
    @boundscheck (m,n) == size(xx) || throw(BoundsError())
    @boundscheck (m,n) == size(yy) || throw(BoundsError())

    if (m+2+ep,n+4) == (muv,nuv)      # u case
        # halo sizes
        xh1 = 1+ep
        xh2 = 1
        yh1 = 2
        yh2 = 2
    elseif (m+4,n+2) == (muv,nuv)    # v case
        # halo sizes 1: left/bottom, 2: right/top
        xh1 = 2
        xh2 = 1         # don't fully understand why this is 1 ...
        yh1 = 1
        yh2 = 1
    else
        throw(BoundsError())
    end

    for j ∈ 1:n
        for i ∈ 1:m
            # floor is not defined for posits...
            xi = Int(floor(Float64(xx[i,j])))   # departure point indices lower left corner within grid cell
            yi = Int(floor(Float64(yy[i,j])))

            k,x0 = clip_or_wrap(xi,i,xx[i,j],muv,xh1,xh2)
            l,y0 = clip(yi,j,yy[i,j],nuv,yh1,yh2)

            uvi[i,j] = bilin(uv[k,l],uv[k+1,l],uv[k,l+1],uv[k+1,l+1],x0,y0)
        end
    end
end

"""Advection of sst/tracer based on the departure points xx,yy via bilinear interpolation.
Departure points are clipped/wrapped to remain within the domain. Boundary conditions either
periodic (wrap around behaviour) or no-flux (no gradient via clipping). Once the respective
4 surrounding grid points are found do bilinear interpolation on the unit square."""
function adv_sst!(ssti::AbstractMatrix,sst::AbstractMatrix,xx::AbstractMatrix,yy::AbstractMatrix)
    m,n = size(ssti)
    @boundscheck (m,n) == size(sst) || throw(BoundsError())
    @boundscheck (m-2*halosstx,n-2*halossty) == size(xx) || throw(BoundsError())
    @boundscheck (m-2*halosstx,n-2*halossty) == size(yy) || throw(BoundsError())

    @inbounds for j ∈ 1:n-2*halossty
        for i ∈ 1:m-2*halosstx

            xi = Int(floor(Float64(xx[i,j])))   # departure point lower left corner
            yi = Int(floor(Float64(yy[i,j])))   # coordinates

            k,x0 = clip_or_wrap(xi,i,xx[i,j],m,halosstx,halosstx)
            l,y0 = clip(yi,j,yy[i,j],n,halossty,halossty)

            ssti[i+halosstx,j+halossty] = bilin(sst[k,l],sst[k+1,l],sst[k,l+1],sst[k+1,l+1],x0,y0)
        end
    end
end

"""Clips the relative lower-left corner index xyi (for both x or y indices) to remain within the domain.
ij is the index of the underlying matrix. xy is the actual coordinate, mn (m or n) the size of the domain,
and h1,h2 are the halo sizes (left/south and right/north)."""
function clip(xyi::Int,ij::Int,xy::Real,mn::Int,h1::Int,h2::Int)

    xyis = xyi+ij+h1                # shifted index (i.e. revert the relative index to an absolute)

    if xyis < 1+h1                  # beyond left/southern boundary
        return 1+h1,zeero           # coordinate is then 0
    elseif xyis > mn-1           # beyond right/northern boundary
        return mn-1,oone         # coordinate is then 1
    else                            # normal case.
        return xyis,xy-xyi          # relative coordinate ∈ [0,1]
    end
end

"""Clips the relative lower-left corner index xyi (for both x or y indices) to remain within the domain.
ij is the index of the underlying matrix. xy is the actual coordinate, mn (m or n) the size of the domain,
and h is the halo size."""
function wrap(xyi::Int,ij::Int,xy::Real,mn::Int,h1::Int,h2::Int)

    xyis = xyi+ij+h1                 # shifted index (i.e. revert the relative index to an absolute)
    xy0 = xy-xyi                     # relative coordinate ∈ [0,1]

    if xyis < 1+h1                   # beyond left/southern boundary
        return xyis+mn-h1-h2,xy0     # coordinate is wrapped around
    elseif xyis > mn-1              # beyond right/northern boundary
        return xyis-mn+h1+h2,xy0     # coordinate is wrapped around
    else                             # normal case.
        return xyis,xy0              # just pass
    end
end

"""Bilinear interpolation on (x,y) in the unit square [0,1]x[0,1].
The values at the corners are f00 = f(0,0), f01 = f(0,1), etc."""
function bilin(f00::T,f10::T,f01::T,f11::T,x::T,y::T) where {T<:AbstractFloat}
    oone = one(T)
    return f00*(oone-x)*(oone-y) + f10*x*(oone-y) + f01*(oone-x)*y + f11*x*y
end

"""Tracer relaxation."""
function tracer_relax!(sst::AbstractMatrix,sst_ref::AbstractMatrix,SSTγ::AbstractMatrix)
    m,n = size(sst)
    @boundscheck (m-2*halosstx,n-2*halossty) == size(sst_ref) || throw(BoundsError())
    @boundscheck (m-2*halosstx,n-2*halossty) == size(SSTγ) || throw(BoundsError())

    @inbounds for j ∈ 1+halossty:n-halossty
        for i ∈ 1+halosstx:m-halosstx
            sst[i,j] += SSTγ[i-halosstx,j-halossty]*(sst_ref[i-halosstx,j-halossty] - sst[i,j])
        end
    end
end

"""Tracer consumption."""
function tracer_consumption!(sst::AbstractMatrix)
    m,n = size(sst)

    @inbounds for j ∈ 1+halossty:n-halossty
        for i ∈ 1+halosstx:m-halosstx
            sst[i,j] += SST_J*(SST0 - sst[i,j])
        end
    end
end

"""Spatially dependent relaxation time scale."""
function sst_γ(x::AbstractVector,y::AbstractVector)
    xx,yy = meshgrid(x,y)

    # convert from days to one over 1/s and include adv time step
    γ0 = dtadvint/(SST_γ0*3600*24)

    x10E = 10*m_per_lat()   # assume Equator: lat/lon equivalence
    γ = γ0/2 .* (1 .- tanh.((xx.-SST_λ0)./SST_λs))
    γ[xx .> x10E] .= 0.0
    return Numtype.(γ)
end

if bc_x == "periodic"
    clip_or_wrap = wrap
else
    clip_or_wrap = clip
end

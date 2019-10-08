""" Extends the matrices u,v,η with a halo of ghost points for boundary conditions. """
function add_halo(u,v,η,sst)
    if VERSION == v"0.6.1"  # cat function was different back then
        u = cat(1,zeros(Numtype,halo,nuy),u,zeros(Numtype,halo,nuy))
        u = cat(2,zeros(Numtype,nux+2*halo,halo),u,zeros(Numtype,nux+2*halo,halo))

        v = cat(1,zeros(Numtype,halo,nvy),v,zeros(Numtype,halo,nvy))
        v = cat(2,zeros(Numtype,nvx+2*halo,halo),v,zeros(Numtype,nvx+2*halo,halo))

        # and haloη for η
        η = cat(1,zeros(Numtype,haloη,ny),η,zeros(Numtype,haloη,ny))
        η = cat(2,zeros(Numtype,nx+2*haloη,haloη),η,zeros(Numtype,nx+2*haloη,haloη))

        sst = cat(1,zeros(Numtype,halosstx,ny),sst,zeros(Numtype,halosstx,ny))
        sst = cat(2,zeros(Numtype,nx+2*halosstx,halossty),sst,zeros(Numtype,nx+2*halosstx,halossty))

    else
        # use halo = number of halo rows/columns to either side for u,v
        u = cat(zeros(Numtype,halo,nuy),u,zeros(Numtype,halo,nuy),dims=1)
        u = cat(zeros(Numtype,nux+2*halo,halo),u,zeros(Numtype,nux+2*halo,halo),dims=2)

        v = cat(zeros(Numtype,halo,nvy),v,zeros(Numtype,halo,nvy),dims=1)
        v = cat(zeros(Numtype,nvx+2*halo,halo),v,zeros(Numtype,nvx+2*halo,halo),dims=2)

        # and haloη for η
        η = cat(zeros(Numtype,haloη,ny),η,zeros(Numtype,haloη,ny),dims=1)
        η = cat(zeros(Numtype,nx+2*haloη,haloη),η,zeros(Numtype,nx+2*haloη,haloη),dims=2)

        sst = cat(zeros(Numtype,halosstx,ny),sst,zeros(Numtype,halosstx,ny),dims=1)
        sst = cat(zeros(Numtype,nx+2*halosstx,halossty),sst,zeros(Numtype,nx+2*halosstx,halossty),dims=2)
    end

    ghost_points!(u,v,η)
    ghost_points_sst!(sst)
    return u,v,η,sst
end

""" Copy ghost points for u from inside to the halo in the nonperiodic case. """
function ghost_points_u_nonperiodic!(u)

    # kinematic boundary condition
    # these values shouldn't change during a time step
    # u[1,:] = zeero
    # u[2,:] = zeero
    # u[end-1,:] = zeero
    # u[end,:] = zeero

    # tangential boundary condition
    @views @inbounds u[:,1] .= one_minus_α*u[:,4]
    @views @inbounds u[:,2] .= one_minus_α*u[:,3]
    @views @inbounds u[:,end-1] .= one_minus_α*u[:,end-2]
    @views @inbounds u[:,end] .= one_minus_α*u[:,end-3]
end

""" Copy ghost points for u from inside to the halo in the periodic case. """
function ghost_points_u_periodic!(u)

    # periodic bc
    @views @inbounds u[1,:] .= u[end-3,:]
    @views @inbounds u[2,:] .= u[end-2,:]
    @views @inbounds u[end-1,:] .= u[3,:]
    @views @inbounds u[end,:] .= u[4,:]

    # tangential bc
    @views @inbounds u[:,1] .= one_minus_α*u[:,4]
    @views @inbounds u[:,2] .= one_minus_α*u[:,3]
    @views @inbounds u[:,end-1] .= one_minus_α*u[:,end-2]
    @views @inbounds u[:,end] .= one_minus_α*u[:,end-3]
end

""" Copy ghost points for v from inside to the halo in the nonperiodic case. """
function ghost_points_v_nonperiodic!(v)

    # kinematic boundary condition
    # these values shouldn't change during a time step
    # v[:,1] .= zeero
    # v[:,2] .= zeero
    # v[:,end-1] .= zeero
    # v[:,end] .= zeero

    # tangential boundary condition
    @views @inbounds v[1,:] .= one_minus_α*v[4,:]
    @views @inbounds v[2,:] .= one_minus_α*v[3,:]
    @views @inbounds v[end-1,:] .= one_minus_α*v[end-2,:]
    @views @inbounds v[end,:] .= one_minus_α*v[end-3,:]
end

""" Copy ghost points for v from inside to the halo in the periodic case. """
function ghost_points_v_periodic!(v)

    # kinematic boundary condition
    # these values shouldn't change during a time step
    # v[:,1] .= zeero
    # v[:,2] .= zeero
    # v[:,end-1] .= zeero
    # v[:,end] .= zeero

    # tangential boundary condition
    @views @inbounds v[1,:] .= v[end-3,:]
    @views @inbounds v[2,:] .= v[end-2,:]
    @views @inbounds v[end-1,:] .= v[3,:]
    @views @inbounds v[end,:] .= v[4,:]
end

""" Copy ghost points for η from inside to the halo in the nonperiodic case. """
function ghost_points_η_nonperiodic!(η)

    # assume no gradients of η across solid boundaries
    # the 4 corner points are copied twice, but it's faster!
    @views @inbounds η[1,:] .= η[2,:]
    @views @inbounds η[end,:] .= η[end-1,:]

    @views @inbounds η[:,1] .= η[:,2]
    @views @inbounds η[:,end] .= η[:,end-1]
end

""" Copy ghost points for η from inside to the halo in the periodic case. """
function ghost_points_η_periodic!(η)

    # corner points are copied twice, but it's faster!
    @views @inbounds η[1,:] .= η[end-1,:]
    @views @inbounds η[end,:] .= η[2,:]

    @views @inbounds η[:,1] .= η[:,2]
    @views @inbounds η[:,end] .= η[:,end-1]
end

""" Copy ghost points for η from inside to the halo in the nonperiodic case. """
function ghost_points_sst_nonperiodic!(sst)

    # assume no gradients of η across solid boundaries
    # the 4 corner points are copied twice, but it's faster!
    for i ∈ 1:halosstx
        @views @inbounds sst[i,:] .= sst[halosstx+1,:]
        @views @inbounds sst[end-i+1,:] .= sst[end-halosstx,:]
    end

    for j ∈ 1:halossty
        @views @inbounds sst[:,j] .= sst[:,halossty+1]
        @views @inbounds sst[:,end-j+1] .= sst[:,end-halossty]
    end
end

""" Copy ghost points for η from inside to the halo in the periodic case. """
function ghost_points_sst_periodic!(sst)

    # corner points are copied twice, but it's faster!
    for i ∈ 1:halosstx
        @views @inbounds sst[i,:] .= sst[end-2*halosstx+i,:]
        @views @inbounds sst[end-halosstx+i,:] .= sst[halosstx+i,:]
    end

    for j ∈ 1:halossty
        @views @inbounds sst[:,j] .= sst[:,halossty+1]
        @views @inbounds sst[:,end-j+1] .= sst[:,end-halossty]
    end
end

# gather and rename functions for convenience
function ghost_points_periodic!(u,v,η)

    ghost_points_u_periodic!(u)
    ghost_points_v_periodic!(v)
    ghost_points_η_periodic!(η)

end

function ghost_points_nonperiodic!(u,v,η)

    ghost_points_u_nonperiodic!(u)
    ghost_points_v_nonperiodic!(v)
    ghost_points_η_nonperiodic!(η)

end

# pick the right set of functions depending on boundary conditions
if bc_x == "periodic"
    ghost_points! = ghost_points_periodic!
    ghost_points_sst! = ghost_points_sst_periodic!
else
    ghost_points! = ghost_points_nonperiodic!
    ghost_points_sst! = ghost_points_sst_nonperiodic!
end

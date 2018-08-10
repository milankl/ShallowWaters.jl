function add_halo(u,v,η)

    u = cat(1,zeros(Numtype,halo,nuy),u,zeros(Numtype,halo,nuy))
    u = cat(2,zeros(Numtype,nux+2*halo,halo),u,zeros(Numtype,nux+2*halo,halo))

    v = cat(1,zeros(Numtype,halo,nvy),v,zeros(Numtype,halo,nvy))
    v = cat(2,zeros(Numtype,nvx+2*halo,halo),v,zeros(Numtype,nvx+2*halo,halo))

    # halo is always 1 for η
    η = cat(1,zeros(Numtype,1,ny),η,zeros(Numtype,1,ny))
    η = cat(2,zeros(Numtype,nx+2,1),η,zeros(Numtype,nx+2,1))

    ghost_points!(u,v,η)
    return u,v,η
end

function ghost_points_u_nonperiodic!(u)

    # kinematic boundary condition
    # these values shouldn't change during a time step
    # u[1,:] = 0.
    # u[2,:] = 0.
    # u[end-1,:] = 0.
    # u[end,:] = 0.

    @views u[:,1] .= one_minus_α*u[:,4]
    @views u[:,2] .= one_minus_α*u[:,3]
    @views u[:,end-1] .= one_minus_α*u[:,end-2]
    @views u[:,end] .= one_minus_α*u[:,end-3]
end

function ghost_points_u_periodic!(u)

    # periodic bc
    @views u[1,:] .= u[end-3,:]
    @views u[2,:] .= u[end-2,:]
    @views u[end-1,:] .= u[3,:]
    @views u[end,:] .= u[4,:]

    # tangential bc
    @views u[:,1] .= one_minus_α*u[:,4]
    @views u[:,2] .= one_minus_α*u[:,3]
    @views u[:,end-1] .= one_minus_α*u[:,end-2]
    @views u[:,end] .= one_minus_α*u[:,end-3]
end

function ghost_points_v_nonperiodic!(v)

    # kinematic boundary condition
    # these values shouldn't change during a time step
    # v[:,1] = zeero
    # v[:,2] = zeero
    # v[:,end-1] = zeero
    # v[:,end] = zeero

    @views v[1,:] .= one_minus_α*v[4,:]
    @views v[2,:] .= one_minus_α*v[3,:]
    @views v[end-1,:] .= one_minus_α*v[end-2,:]
    @views v[end,:] .= one_minus_α*v[end-3,:]
end

function ghost_points_v_periodic!(v)

    # kinematic boundary condition
    # these values shouldn't change during a time step
    # v[:,1] = zeero
    # v[:,2] = zeero
    # v[:,end-1] = zeero
    # v[:,end] = zeero

    @views v[1,:] .= v[end-3,:]
    @views v[2,:] .= v[end-2,:]
    @views v[end-1,:] .= v[3,:]
    @views v[end,:] .= v[4,:]
end

function ghost_points_η_nonperiodic!(η)

    # corner points are copied twice, but it's faster!
    @views η[1,:] .= η[2,:]
    @views η[end,:] .= η[end-1,:]

    @views η[:,1] .= η[:,2]
    @views η[:,end] .= η[:,end-1]
end

function ghost_points_η_periodic!(η)

    # corner points are copied twice, but it's faster!
    @views η[1,:] .= η[end-1,:]
    @views η[end,:] .= η[2,:]

    @views η[:,1] .= η[:,2]
    @views η[:,end] .= η[:,end-1]
end

## GATHER AND RENAME FUNCTIONS FOR CONVENIENCE

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

# pick the right set of functions depending on boundary conditions and halo size

if bc_x == "periodic"
    ghost_points! = ghost_points_periodic!
else
    ghost_points! = ghost_points_nonperiodic!
end

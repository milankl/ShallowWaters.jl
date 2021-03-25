@with_kw struct Grid{T<:AbstractFloat,Tprog<:AbstractFloat}

    # Parameters taken from Parameter struct
    nx::Int                             # number of grid cells in x-direction
    Lx::Real                            # length of the domain in x-direction [m]
    L_ratio::Real                       # Domain aspect ratio of Lx/Ly
    bc::String                          # boundary condition, "periodic" or "nonperiodic"
    g::Real                             # gravitational acceleration [m/s]
    H::Real                             # layer thickness at rest [m]
    cfl::Real                           # CFL number (1.0 recommended for RK4, 0.6 for RK3)
    Ndays::Real                         # number of days to integrate for
    nstep_diff::Int                     # diffusive terms every nstep_diff time steps.
    nstep_advcor::Int                   # nonlinear terms every nstep_advcor time steps.
    Uadv::Real                          # Velocity scale [m/s] for tracer advection
    output_dt::Real                     # output time step [hours]
    ω::Real                             # Earth's angular frequency [s^-1]
    ϕ::Real                             # central latitue of the domain (for coriolis) [°]
    R::Real                             # Earth's radius [m]
    scale::Real                         # multiplicative scale for momentum equations [1]

    # DOMAIN SIZES
    Δ::Real=Lx / nx                         # grid spacing
    ny::Int=Int(round(Lx / L_ratio / Δ))    # number of grid cells in y-direction
    Ly::Real=ny * Δ                         # length of domain in y-direction

    # NUMBER OF GRID POINTS
    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = ny-1                                     # v-grid in y-direction
    nqx::Int = if (bc == "periodic") nx else nx+1 end   # q-grid in x-direction
    nqy::Int = ny+1                                     # q-grid in y-direction

    # TOTAL NUMBER OF GRID POINTS
    nT::Int = nx*ny                     # T-grid
    nu::Int = nux*nuy                   # u-grid
    nv::Int = nvx*nvy                   # v-grid
    nq::Int = nqx*nqy                   # q-grid

    # GRID VECTORS
    x_T::AbstractVector = Δ*Array(1:nx) .- Δ/2
    y_T::AbstractVector = Δ*Array(1:ny) .- Δ/2
    x_u::AbstractVector = if (bc == "periodic") Δ*Array(0:nx-1) else Δ*Array(1:nx-1) end
    y_u::AbstractVector = y_T
    x_v::AbstractVector = x_T
    y_v::AbstractVector = Δ*Array(1:ny-1)
    x_q::AbstractVector = if bc == "periodic" x_u else Δ*Array(1:nx+1) .- Δ end
    y_q::AbstractVector = Δ*Array(1:ny+1) .- Δ

    # HALO SIZES
    halo::Int=2                         # halo size for u,v (Biharmonic stencil requires 2)
    haloη::Int=1                        # halo size for η
    halosstx::Int=1                     # halo size for tracer sst in x
    halossty::Int=0                     # halo size for tracer sst in y

    # EDGE POINT (1 = yes, 0 = no)
    ep::Int = if bc == "periodic" 1 else 0 end  # is there a u-point on the left edge?

    # GRID VECTORS WITH HALO
    x_T_halo::AbstractVector = Δ*Array(0:nx+1) .- Δ/2
    y_T_halo::AbstractVector = Δ*Array(0:ny+1) .- Δ/2
    x_u_halo::AbstractVector = if (bc == "periodic") Δ*Array(-2:nx+1) else Δ*Array(-1:nx+1) end
    y_u_halo::AbstractVector = Δ*Array(-1:ny+2) .- Δ/2
    x_v_halo::AbstractVector = Δ*Array(-1:nx+2) .- Δ/2
    y_v_halo::AbstractVector = Δ*Array(-1:ny+1)
    x_q_halo::AbstractVector = if bc == "periodic" x_u_halo else Δ*Array(-1:nx+3) .- Δ end
    y_q_halo::AbstractVector = Δ*Array(-1:ny+3) .- Δ

    # TIME STEPS
    c::Real = √(g*H)                            # shallow water gravity wave speed
    dtint::Int = Int(floor(cfl*Δ/c))            # dt converted to Int
    nt::Int = Int(ceil(Ndays*3600*24/dtint))    # number of time steps to integrate
    dt::T = T(dtint)                            # time step [s]
    Δt::T = T(dtint/Δ)                          # time step divided by grid spacing [s/m]
    Δt_diff::Tprog = Tprog(nstep_diff*dtint/Δ)  # time step for diffusive terms

    # TIME STEPS FOR ADVECTION
    nadvstep::Int = max(1,Int(floor(Δ/Uadv/dtint)))         # advection each n time steps
    nadvstep_half::Int = nadvstep ÷ 2                       # nadvstep ÷ 2
    dtadvint::Int = nadvstep*dtint                          # advection time step [s]
    # divide by scale here to undo the scaling in u,v for tracer advection
    dtadvu::T = T(dtadvint*nx/Lx/scale)                     # Rescaled advection time step for u [s/m]
    dtadvv::T = T(dtadvint*ny/Ly/scale)                     # Rescaled advection time step for v [s/m]
    half_dtadvu::T = T(dtadvint*nx/Lx/2/scale)              # dtadvu/2
    half_dtadvv::T = T(dtadvint*ny/Ly/2/scale)              # dtadvv/2

    # N TIME STEPS FOR OUTPUT
    nout::Int = max(1,Int(floor(output_dt*3600/dtint)))     # output every n time steps
    nout_total::Int = (nt ÷ nout)+1                         # total number of output time steps
    t_vec::AbstractVector = Array(0:nout_total-1)*dtint     # time vector

    # CORIOLIS
    f₀::Float64 = coriolis_at_lat(ω,ϕ)                      # Coriolis parameter
    β::Float64 = β_at_lat(ω,R,ϕ)                            # Derivate of Coriolis parameter wrt latitude
    # scale only f_q as it's used for non-linear advection
    f_q::Array{T,2} = T.(scale*Δ*(f₀ .+ β*(yy_q(bc,x_q_halo,y_q_halo) .- Ly/2)))  # same on the q-grid
    # f_u, f_v are only used for linear dynamics (scaling implicit)
    f_u::Array{T,2} = T.(Δ*(f₀ .+ β*(meshgrid(x_u,y_u)[2] .- Ly/2)))        # f = f₀ + βy on the u-grid
    f_v::Array{T,2} = T.(Δ*(f₀ .+ β*(meshgrid(x_v,y_v)[2] .- Ly/2)))        # same on the v-grid
end

"""Helper function to create yy_q based on the boundary condition bc."""
function yy_q(bc::String,x_q_halo::AbstractVector,y_q_halo::AbstractVector)
    if bc == "periodic"
        # points on the right edge needed too
        _,yy_q = meshgrid(x_q_halo[3:end-1],y_q_halo[3:end-2])
    else
        _,yy_q = meshgrid(x_q_halo[3:end-2],y_q_halo[3:end-2])
    end

    return yy_q
end


"""Generator function for the Grid struct."""
function Grid{T,Tprog}(P::Parameter) where {T<:AbstractFloat,Tprog<:AbstractFloat}
    @unpack nx,Lx,L_ratio = P
    @unpack bc,g,H,cfl = P
    @unpack Ndays,nstep_diff,nstep_advcor = P
    @unpack Uadv,output_dt = P
    @unpack ϕ,ω,R,scale = P

    return Grid{T,Tprog}(nx=nx,Lx=Lx,L_ratio=L_ratio,bc=bc,g=g,H=H,cfl=cfl,Ndays=Ndays,
                nstep_diff=nstep_diff,nstep_advcor=nstep_advcor,Uadv=Uadv,output_dt=output_dt,
                ϕ=ϕ,ω=ω,R=R,scale=scale)
end

"""Meter per 1 degree of latitude (or longitude at the equator)."""
function m_per_lat(R::Real)
    return 2π*R/360.
end

"""Coriolis parameter f [1/s] at latitude ϕ [°] given Earth's rotation ω [1/s]."""
function coriolis_at_lat(ω::Real,ϕ::Real)
    return 2*ω*sind(ϕ)
end

"""Coriolis parameter's derivative β wrt latitude [(ms)^-1] at latitude ϕ, given
Earth's rotation ω [1/s] and radius R [m]."""
function β_at_lat(ω::Real,R::Real,ϕ::Real)
    return 2*ω/R*cosd(ϕ)
end

"""Similar to the numpy meshgrid function:
repeats x length(y)-times and vice versa. Returns two matrices xx,yy of same shape so that
each row of xx is x and each column of yy is y."""
function meshgrid(x::AbstractVector,y::AbstractVector)
    m,n = length(x),length(y)

    # preallocate preserving the data type of x,y
    xx = zeros(typeof(x[1]),m,n)
    yy = zeros(typeof(y[1]),m,n)

    for i in 1:m
        xx[i,:] .= x[i]
    end

    for i in 1:n
        yy[:,i] .= y[i]
    end

    return xx,yy
end

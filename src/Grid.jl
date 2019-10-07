mutable struct Grid{T<:AbstractFloat}

    # DOMAIN SIZES
    Δ::Real                             # grid spacing
    nx::Int                             # number of grid cells in x-direction
    ny::Int                             # number of grid cells in y-direction
    Lx::Real                            # length of domain in x-direction
    Ly::Real                            # length of domain in y-direction
    L_ratio::Real                       # aspect ratio of x/y

    # NUMBER OF GRID POINTS
    nux::Int                            # u-grid in x-direction
    nuy::Int                            # u-grid in y-direction
    nvx::Int                            # v-grid in x-direction
    nvy::Int                            # v-grid in y-direction
    nqx::Int                            # q-grid in x-direction
    nqy::Int                            # q-grid in y-direction

    # TOTAL NUMBER OF GRID POINTS
    nT::Int                             # T-grid
    nu::Int                             # u-grid
    nv::Int                             # v-grid
    nq::Int                             # q-grid

    # GRID VECTORS
    x_T::AbstractVector
    y_T::AbstractVector
    x_u::AbstractVector
    y_u::AbstractVector
    x_v::AbstractVector
    y_v::AbstractVector
    x_q::AbstractVector
    y_q::AbstractVector

    # HALO SIZES
    halo::Int                           # halo size for u,v (Biharmonic stencil requires 2)
    haloη::Int                          # halo size for η
    halosstx::Int                       # halo size for tracer sst in x
    halossty::Int                       # halo size for tracer sst in y

    # EDGE POINT (1 = yes, 0 = no)
    ep::Int                             # is there a u-point on the left edge?
                                        # determined by boundary condition

    # GRID VECTORS WITH HALO
    x_T_halo::AbstractVector
    y_T_halo::AbstractVector
    x_u_halo::AbstractVector
    y_u_halo::AbstractVector
    x_v_halo::AbstractVector
    y_v_halo::AbstractVector
    x_q_halo::AbstractVector
    y_q_halo::AbstractVector

    # TIME STEPS
    dt::T                               # time step [s]
    Δt::T                               # time step divided by grid spacing [s/m]
    Δt_diff::T                          # time step for diffusive terms
    dtint::Int                          # dt converted to Int
    nt::Int                             # number of time steps to integrate

    # TIME STEPS FOR ADVECTION
    dtadvu::T                           # Rescaled advection time step for u [s/m]
    dtadvv::T                           # Rescaled advection time step for v [s/m]
    dtadvint::Int                       # advection time step [s]
    nadvstep::Int                       # advection each n time steps
    nadvstep_half::Int                  # nadvstep ÷ 2

    # N TIME STEPS FOR OUTPUT
    nout::Int                           # output every n time steps
    nout_total::Int                     # total number of output time steps
    t_vec::AbstractVector               # time vector

    # CORIOLIS
    f₀::Float64                         # Coriolis parameter
    β::Float64                          # Derivate of Coriolis parameter wrt latitude
    f_u::Array{T,2}                     # f = f₀ + βy on the u-grid
    f_v::Array{T,2}                     # same on the v-grid
    f_q::Array{T,2}                     # same on the q-grid

    # INITIALISE WITH ARBITRARY VALUES FROM MEMORY
    Grid{T}() where T = new{T}()
end

""" Creates the Grid struct based on the parameter struct P."""
function Grid{T}(P::Parameter) where {T<:AbstractFloat}

    @unpack nx,Lx,L_ratio = P
    @unpack bc = P
    @unpack g,H = P
    @unpack cfl = P
    @unpack Ndays,nstep_diff = P
    @unpack Uadv = P

    # INITIALISE
    G = Grid{T}()
    G.nx,G.Lx,G.L_ratio = nx,Lx,L_ratio     # copy values from P to G

    # DOMAIN SIZE - based on the desired aspect ratio L_ratio
    Δ = Lx / nx                             # grid spacing, same in x,y
    ny = Int(round(Lx / L_ratio / Δ))       # number of grid cells in y-direction
    Ly = ny * Δ                             # length of domain in y-direction
    G.Δ,G.ny,G.Ly = Δ,ny,Ly                 # pack

    # NUMBER OF GRID POINTS
    G.nux = if (bc == "periodic") nx else nx-1 end
    G.nuy = ny
    G.nvx,G.nvy = nx,ny-1
    G.nqx = if (bc == "periodic") nx else nx+1 end
    G.nqy = ny+1

    # TOTAL NUMBER of T,u,v,q-points
    G.nT = nx*ny
    G.nu = G.nux*G.nuy
    G.nv = G.nvx*G.nvy
    G.nq = G.nqx*G.nqy

    # GRID VECTORS
    G.x_T = Δ*Array(1:nx) .- Δ/2
    G.y_T = Δ*Array(1:ny) .- Δ/2

    G.x_u = if (bc == "periodic") Δ*Array(0:nx-1) else Δ*Array(1:nx-1) end
    G.y_u = G.y_T

    G.x_v = G.x_T
    G.y_v = Δ*Array(1:ny-1)

    G.x_q = if bc == "periodic" G.x_u else Δ*Array(1:nx+1) .- Δ end
    G.y_q = Δ*Array(1:ny+1) .- Δ

    # EDGE POINT - used in some functions of rhs.jl to avoid an if
    G.ep = if bc == "periodic" 1 else 0 end

    # HALO SIZE
    G.halo = 2          # halo size for u,v (Biharmonic stencil requires 2)
    G.haloη = 1         # halo size for η
    G.halosstx = 1      # halo size for tracer sst in x
    G.halossty = 0      # halo size for tracer sst in y

    # GRID VECTORS WITH HALO
    G.x_T_halo = Δ*Array(0:nx+1) .- Δ/2
    G.y_T_halo = Δ*Array(0:ny+1) .- Δ/2

    G.x_u_halo = if (bc == "periodic") Δ*Array(-2:nx+1) else Δ*Array(-1:nx+1) end
    G.y_u_halo = Δ*Array(-1:ny+2) .- Δ/2

    G.x_v_halo = Δ*Array(-1:nx+2) .- Δ/2
    G.y_v_halo = Δ*Array(-1:ny+1)

    G.x_q_halo = if bc == "periodic" G.x_u_halo else Δ*Array(-1:nx+3) .- Δ end
    G.y_q_halo = Δ*Array(-1:ny+3) .- Δ

    # TIME STEPS - based on CFL to resolve gravity waves
    c = √(g*H)                              # shallow water gravity wave phase speed
    dtint = Int(floor(cfl*Δ/c))             # make the timestep slightly shorter due to floor
    G.nt = Int(ceil(Ndays*3600*24/dtint))   # number of time steps to integrate
    G.dt = T(dtint)                         # convert to number format T
    G.Δt = T(dtint/Δ)                       # timestep combined with grid spacing Δ
    G.Δt_diff = T(nstep_diff*dtint/Δ)       # timestep for diffusive terms
    G.dtint = dtint                         # pack

    G.nout = Int(floor(P.output_dt*3600/dtint)) # output every nout time steps
    G.nout_total = (G.nt ÷ G.nout)+1            # total number of output time steps
    G.t_vec = Array(0:G.nout_total-1)*dtint     # time vector for output

    # ADVECTION TIME STEPS
    # round down to make sure nadvstep is integer, at least 1
    nadvstep = max(1,Int(floor(Δ/Uadv/dtint)))  # advection every n time steps
    G.nadvstep_half = nadvstep ÷ 2              # for mid-point calculation
    dtadvint = nadvstep*dtint                   # time step for advection
    G.nadvstep = nadvstep
    G.dtadvint = dtadvint                       # pack

    # For dimensionless grid: rescaled advection time step [s/m]
    G.dtadvu = T(dtadvint*nx/Lx)
    G.dtadvv = T(dtadvint*ny/Ly)

    # CORIOLIS
    @unpack ω,ϕ,R = P
    f₀ = coriolis_at_lat(ω,ϕ)
    β = β_at_lat(ω,R,ϕ)

    xx_u,yy_u = meshgrid(G.x_u,G.y_u)
    xx_v,yy_v = meshgrid(G.x_v,G.y_v)

    if bc == "periodic"
        # points on the right edge needed too
        xx_q,yy_q = meshgrid(G.x_q_halo[3:end-1],G.y_q_halo[3:end-2])
    else
        xx_q,yy_q = meshgrid(G.x_q,G.y_q)
    end

    G.f_u = T.(Δ*(f₀ .+ β*(yy_u .- Ly/2)))
    G.f_v = T.(Δ*(f₀ .+ β*(yy_v .- Ly/2)))
    G.f_q = T.(Δ*(f₀ .+ β*(yy_q .- Ly/2)))

    G.f₀ = f₀
    G.β = β

    return G
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

"""Based on the desired aspect ratio L_ratio of the domain returns the number of
grid cells in y-direction and the length of that side of the domain."""
function domain_ratio(nx::Int,Lx::Real,L_ratio::Real)
    Δ = Lx / nx                        # grid spacing
    ny = Int(round(Lx / L_ratio / Δ))  # number of grid cells in y-direction
    Ly = ny * Δ                        # length of domain in y-direction
    return Δ,ny,Ly
end

"""Similar to the numpy meshgrid function:
repeats x length(y)-times and vice versa. Returns two matrices xx,yy of same shape so that
each row of xx is x and each column of yy is y."""
function meshgrid(x,y)
    # preallocate preserving the data type of x,y
    xx = zeros(typeof(x[1]),length(x),length(y))
    yy = zeros(typeof(y[1]),length(x),length(y))

    for i in 1:length(x)
        xx[i,:] .= x[i]
    end

    for i in 1:length(y)
        yy[:,i] .= y[i]
    end

    return xx,yy
end

""" Returns the time step dt,Δt,dtint and the total number of steps to integrate nt
based on the CFL criterion to resolve shallow water gravity waves. """
function timestep()
    # shallow water gravity wave phase speed
    c_phase = sqrt(gravity*water_depth)

    # the model timestep dt based on cfl stability criterion to resolve gravity waves
    # converting to integer, i.e. rounding up (ceil)
    dtint = Int(floor(cfl*Δ/c_phase))      # make the timestep slightly shorter due to floor
    nt = Int(ceil(Ndays*3600*24/dtint))    # number of time steps to integrate
    dt = Numtype(dtint)                    # convert to Numtype for multiplication of the RHS
    Δt = Numtype(dtint/Δ)                  # timestep combined with grid spacing Δ
    return dt,Δt,dtint,nt
end

""" Returns the tracer advection time step dtadv and the number of timesteps nadvstep after which
one evaluation of the tracer advection is computed. """
function adv_timestep()
    # round down to make sure nadvstep is integer, at least 1
    nadvstep = max(1,Int(floor(Δ/Uadv/dtint)))
    nadvstep_half = nadvstep ÷ 2    #
    # recompute the tracer advection time step to fit the rounding
    dtadvint = nadvstep*dtint
    dtadvu = Numtype(dtadvint*nx/Lx)    # [s/m] for dimensionless advection grid
    dtadvv = Numtype(dtadvint*ny/Ly)
    return dtadvu,dtadvv,dtadvint,nadvstep,nadvstep_half
end

const Δ,ny,Ly = domain_ratio(nx,Lx,L_ratio)

# number of grid points for u,v,q-grid in either x or y-direction
const nux = if (bc_x == "periodic") nx else nx-1 end
const nuy = ny
const nvx,nvy = nx,ny-1
const nqx = if (bc_x == "periodic") nx else nx+1 end
const nqy = ny+1

# total number of T,u,v,q-points
const nT = nx*ny
const nu = nux*nuy
const nv = nvx*nvy
const nq = nqx*nqy

# grid vectors
const x_T = Δ*Array(1:nx) .- Δ/2
const y_T = Δ*Array(1:ny) .- Δ/2

const x_u = if (bc_x == "periodic") Δ*Array(0:nx-1) else Δ*Array(1:nx-1) end
const y_u = y_T

const x_v = x_T
const y_v = Δ*Array(1:ny-1)

const x_q = if bc_x == "periodic" x_u else Δ*Array(1:nx+1) .- Δ end
const y_q = Δ*Array(1:ny+1) .- Δ

# halo of ghost points (because of the biharmonic operator) - don't change.
const halo = 2
const haloη = 1
const halosstx = 2
const halossty = 1

# is there a point on the left edge? ep - egde points
# used in some functions of rhs.jl to avoid an if
const ep = if bc_x == "periodic" 1 else 0 end
#TODO possibly different edge point variables needed on subdomains

# halo versions with additional ghost points
const x_T_halo = Δ*Array(0:nx+1) .- Δ/2
const y_T_halo = Δ*Array(0:ny+1) .- Δ/2

const x_u_halo = if (bc_x == "periodic") Δ*Array(-2:nx+1) else Δ*Array(-1:nx+1) end
const y_u_halo = Δ*Array(-1:ny+2) .- Δ/2

const x_v_halo = Δ*Array(-1:nx+2) .- Δ/2
const y_v_halo = Δ*Array(-1:ny+1)

# also two halo for q
const x_q_halo = if bc_x == "periodic" x_u_halo else Δ*Array(-1:nx+3) .- Δ end
const y_q_halo = Δ*Array(-1:ny+3) .- Δ

# matrices of x and y positions with halo (dimensionless - actually indices!)
const xxT,yyT = meshgrid(Numtype.(Array(1:nx)),Numtype.(Array(1:ny)))

# time and output
const dt,Δt,dtint,nt = timestep()
const nout = Int(floor(output_dt*3600/dtint))   # output every nout time steps
const nout_total = (nt ÷ nout)+1                # total number of time steps for output
const t_vec = Array(0:nout_total-1)*dtint       # time vector for output

# advection time step
const dtadvu,dtadvv,dtadvint,nadvstep,nadvstep_half = adv_timestep()
println((nout,nadvstep))

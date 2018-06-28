function domain_ratio(nx::Int,Lx::Real,L_ratio::Real)
    #= Based on the desired aspect ratio L_ratio of the domain returns the number of
    grid cells in y-direction and the length of that side of the domain. =#

    dx = Lx / nx                        # grid spacing
    ny = Int(round(Lx / L_ratio / dx))  # number of grid cells in y-direction
    Ly = ny * dx                        # length of domain in y-direction
    return dx,ny,Ly
end

function meshgrid(x,y)
    #= Similar to the numpy meshgrid function:
    repeats x length(y)-times and vice versa. Returns two matrices xx,yy of same shape so that
    each row of xx is x and each column of yy is y. =#

    # preallocate preserving the data type of x,y
    xx = zeros(typeof(x[1]),length(x),length(y))
    yy = zeros(typeof(y[1]),length(x),length(y))

    for i in 1:length(x)
        xx[i,:] = x[i]
    end

    for i in 1:length(y)
        yy[:,i] = y[i]
    end

    return xx,yy
end

const dx,ny,Ly = domain_ratio(nx,Lx,L_ratio)

# number of grid points for u,v,q-grid in either x or y-direction
const nux = if (bc_x == "periodic") nx else nx-1 end
const nuy = ny
const nvx,nvy = nx,ny-1
const nqx = if (bc_x == "periodic") nx else nx+1 end
const nqy = ny+1

# total number of T,u,v,q-points
const NT = nx*ny
const Nu = nux*nuy
const Nv = nvx*nvy
const Nq = nqx*nqy

# grid vectors
const x_T = Array(1:nx)*dx - dx/2
const y_T = Array(1:ny)*dx - dx/2

const x_u = if (bc_x == "periodic") Array(0:nx-1)*dx else const x_u = Array(1:nx-1)*dx end
const y_u = y_T

const x_v = x_T
const y_v = Array(1:ny-1)*dx

const x_q = if bc_x == "periodic" x_u else Array(1:nx+1)*dx - dx end
const y_q = Array(1:ny+1)*dx - dx

function grid(nx::Int,Lx::AbstractFloat,L_ratio::AbstractFloat)
    dx = Lx / nx
    ny = Int(round(Lx / L_ratio / dx))
    Ly = ny * dx
    return dx,ny,Ly
end

nx = 3
Lx = 3.
L_ratio = 1.
const α = 2.

#TODO depends on the boundary conditions! periodic vs non-periodic
const dx,ny,Ly = grid(nx,Lx,L_ratio)

const NT = nx*ny
const Nu = (nx-1)*ny   # number of u-points
const Nv = nx*(ny-1)   # number of v-points
const Nq = (nx+1)*(ny+1) # number of q-points

const x_T = Array((1:nx-1)*dx - dx/2)
const y_T = Array((1:ny-1)*dx - dx/2)

const y_u = y_T

const x_v = x_T
const y_v = Array((1:ny-1)*dx)
const y_q = Array((1:ny+1)*dx - dx)

if bc_x == "periodic"
    const x_u = Array((1:nx)*dx - dx)
    const x_q = x_u
else
    const x_u = Array((1:nx-1)*dx)
    const x_q = Array((1:nx+1)*dx - dx)
end

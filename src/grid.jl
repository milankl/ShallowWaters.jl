function grid(nx::Int,Lx::Float,L_ratio::Float)
    dx = Lx / nx
    ny = Int(round(Lx / L_ratio / dx))
    Ly = ny * dx
    return dx,ny,Ly
end


#TODO depends on the boundary conditions! periodic vs non-periodic
const NT = nx*ny
const Nu = (nx-1)*ny   # number of u-points
const Nv = nx*(ny-1)   # number of v-points
const Nq = (nx+1)*(ny+1) # number of q-points

const dx,ny,Ly = grid(nx,Lx,L_ratio)
const one_over_dx = 1/dx


#=
    # grid vectors for T-points
    param['x_T'] = np.arange(param['dx']/2.,param['Lx'],param['dx'])
    param['y_T'] = np.arange(param['dy']/2.,param['Ly'],param['dy'])

    # grid vectors for u-points
    param['x_u'] = param['x_T'][:-1] + param['dx']/2.
    param['y_u'] = param['y_T']

    #grid vectors for v-points
    param['x_v'] = param['x_T']
    param['y_v'] = param['y_T'][:-1] + param['dy']/2.

    # grid vectors for q-points
    param['x_q'] = np.arange(0,param['Lx']+param['dx']/2.,param['dx'])
    param['y_q'] = np.arange(0,param['Ly']+param['dy']/2.,param['dy'])
=#

# define constants
const Numtype = Float32

const nx = 64                   # number of grid cells in x-direction
const Lx = 3840e3               # length of the domain in x-direction
const L_ratio = 1               # Domain aspect ratio of Lx/Ly

const g = 10.                   # gravitational acceleration
const H = 500.                  # layer thickness at rest
const ρ = 1e3                   # density

const cfl = 0.9                 # CFL number
const Ndays = 500                # number of days to integrate for

# boundary condtions
const bc_x = "nonperiodic"      # "periodic" or anything else for nonperiodic
const α = 0.                    # lateral boundary condition parameter
                                # 0 free-slip, 0<α<2 partial-slip, 2 no-slip

const run_id = 0                #TODO make automatic
const c_D = 1e-5

const output = 1

const ϕ = 30.                   # central latitue of the domain

# define constants

const nx = 64                   # number of grid cells in x-direction
const Lx = 512e3                # length of the domain in x-direction
const L_ratio = 3               # Domain aspect ratio of Lx/Ly

const g = 10.                   # gravitational acceleration
const H = 500.                  # layer thickness at rest
const ρ = 1e3                   # density

const cfl = 0.9                 # CFL number
const Ndays = 10                # number of days to integrate for

# boundary condtions
const bc_x = "non-periodic"     # 1 for periodic, 0 for non-periodic
const α = 0.                    # lateral boundary condition parameter

const run_id = 0                #TODO make automatic
const c_D = 1e-5

const output = 1

const ϕ = 30.                   # central latitue of the domain
const Numtype = Float32

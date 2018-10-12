# define constants
const Numtype = Float32
#const Numtype = Posit{16,0}

const nx = 200                  # number of grid cells in x-direction
const Lx = 2000e3               # length of the domain in x-direction
const L_ratio = 2               # Domain aspect ratio of Lx/Ly

const gravity = 10.             # gravitational acceleration
const water_depth = 500.        # layer thickness at rest
const ρ = 1e3                   # density

const wind_forcing = "channel"  # "channel", "double_gyre", "shear" or "none"
const Fx0 = 0.12                  # wind stress strength [Pa], default 0.12

const topography_feature = "ridge"
const topofeat_height = 50.      # height of seamount
const topofeat_width = 100e3    # horizontal scale [m] of the seamount

const surface_forcing = false   # or true
const t_relax = 5.              # time scale of the interface_relaxation [days]
const η_refh = 5.               # height difference [m] of the interface relaxation profile
const η_refw = 50e3             # width [m] of the tangent used for the interface relaxation

const cfl = 0.9                 # CFL number
const Ndays = 400               # number of days to integrate for

# boundary condtions
const bc_x = "periodic"         # "periodic" or anything else for nonperiodic
const lbc = 2.                  # lateral boundary condition parameter
                                # 0 free-slip, 0<lbc<2 partial-slip, 2 no-slip

const adv_scheme = "ArakawaHsu" # "Sadourny" or "ArakawaHsu"

const drag = 1e-5               # bottom drag coefficient [dimensionless]
const c_smag = 0.15             # Smagorinsky coefficient [dimensionless]

const output = 1                # 1 for nc output 0 for none
const output_dt = 6             # output time step in hours
const outpath = "/network/aopp/cirrus/pred/kloewer/julsdata/"

const initial_cond = "rest"   # "rest" or "ncfile"
const init_run_id = 0           # only for starting from ncfile

const ϕ = 30.                   # central latitue of the domain

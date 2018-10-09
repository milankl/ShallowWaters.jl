# define constants
const Numtype = Float32
#const Numtype = Posit{16,0}

const nx = 200                  # number of grid cells in x-direction
const Lx = 2000e3               # length of the domain in x-direction
const L_ratio = 2               # Domain aspect ratio of Lx/Ly

const gravity = 10.             # gravitational acceleration
const water_depth = 500.        # layer thickness at rest
const ρ = 1e3                   # density

const wind_forcing = "double_gyre"
const Fx0 = 0.12                  # wind stress strength [Pa], default 0.12

const topography_feature = "ridge"
const seamount_height = 50.      # height of seamount
const seamount_width = 100e3    # horizontal scale [m] of the seamount

const t_relax = 5.              # time scale of the interface_relaxation [days]
const η_refh = 5.               # height difference [m] of the interface relaxation profile
const η_refw = 50e3            # width [m] of the tangent used for the interface relaxation

const cfl = 0.9                 # CFL number
const Ndays = 200               # number of days to integrate for

# boundary condtions
const bc_x = "periodic"      # "periodic" or anything else for nonperiodic
const lbc = 0.                  # lateral boundary condition parameter
                                # 0 free-slip, 0<lbc<2 partial-slip, 2 no-slip

const adv_scheme = "ArakawaHsu" # "Sadourny" or "ArakawaHsu"

const drag = 1e-5               # bottom drag coefficient
const c_smag = 0.15             # Smagorinsky coefficient

const output = 1                # 1 for nc output 0 for none
const output_dt = 3             # output time step in hours
const outpath = "/local/kloewer/julsdata/"

const initial_cond = "rest"     # "rest" or "ncfile"
const init_run_id = 2           # only for starting from ncfile

const ϕ = 30.                   # central latitue of the domain

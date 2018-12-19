# NUMBER FORMAT OPTIONS
const Numtype = Float32
#const Numtype = Posit{16,2}
#const Numtype = Main.FiniteFloats.Finite16
#const Numtype = BigFloat
#setprecision(7)

# DOMAIN RESOLUTION AND RATIO
const nx = 100                  # number of grid cells in x-direction
const Lx = 2000e3               # length of the domain in x-direction
const L_ratio = 2               # Domain aspect ratio of Lx/Ly

# PHYSICAL CONSTANTS
const gravity = 10.             # gravitational acceleration
const water_depth = 500.        # layer thickness at rest
const ρ = 1e3                   # density
const ϕ = 30.                    # central latitue of the domain (for coriolis)
const ω = 2π/(24*3600)          # Earth's angular frequency [s^-1]
const R = 6.371e6               # Earth's radius [m]

# WIND FORCING OPTIONS
const wind_forcing = "channel"  # "channel", "double_gyre", "shear" or "none"
const Fx0 = 0.12                # wind stress strength [Pa], default 0.12

# BOTTOM TOPOGRAPHY OPTIONS
const topography_feature = "ridges" # "ridge", "seamount", "flat"
const topofeat_height = 80.      # height of seamount
const topofeat_width = 200e3    # horizontal scale [m] of the seamount

# NEWTONIAN COOLING OPTIONS
const surface_forcing = false   # or true
const t_relax = 100.              # time scale of the interface_relaxation [days]
const η_refh = 5.               # height difference [m] of the interface relaxation profile
const η_refw = 100e3             # width [m] of the tangent used for the interface relaxation

# TIME STEPPING OPTIONS
const RKo = 4                   # Order of the RK time stepping scheme (3 or 4)
const cfl = 1.0                 # CFL number (1.0 recommended for RK4, 0.6 for RK3)
const Ndays = 100            # number of days to integrate for

# BOUNDARY CONDITION OPTIONS
const bc_x = "periodic"         # "periodic" or anything else for nonperiodic
const lbc = 1.                  # lateral boundary condition parameter
                                # 0 free-slip, 0<lbc<2 partial-slip, 2 no-slip

# MOMENTUM ADVECTION OPTIONS
const adv_scheme = "ArakawaHsu"   # "Sadourny" or "ArakawaHsu"

# BOTTOM FRICTION OPTIONS
const bottom_friction = "linear" # "linear" or "quadratic"
const drag = 1e-5               # bottom drag coefficient [dimensionless] for quadratic
const τdrag = 300.               # bottom drag coefficient [days] for linear

# DIFFUSION OPTIONS
const diffusion = "Constant"    # "Smagorinsky" or "Constant", biharmonic in both cases
const ν_const = 500             # [m^2/s] scaling constant for Constant biharmonic diffusion
const c_smag = 0.15             # Smagorinsky coefficient [dimensionless]

# TRACER ADVECTION
const tracer_advection = true   # "true" or "false"
const tracer_relaxation = true  # "true" or "false"
const injection_region = "rect"  # "west", "south" or "rect"
const sstrestart = false        # start from previous sst file
const Uadv = 0.5               # Velocity scale [m/s] for tracer advection
const SSTmax = 30.               # tracer (sea surface temperature) max for restoring
const SSTmin = 0.               # tracer (sea surface temperature) min for restoring
const τSST = 600.                # tracer restoring time scale [days]
const SSTw = 1000e3               # width [m] of the tangent used for the IC and interface relaxation
const SSTϕ = 0.5                # latitude/longitude ∈ [0,1] of sst edge

# OUTPUT OPTIONS
const output = true                # 1 for nc output 0 for none
const output_vars = ["u","v","eta","sst"]
const output_dt = 6             # output time step in hours
#const outpath = "/network/aopp/chaos/pred/kloewer/julsdata/ssthr/"
const outpath = "/Users/milan/phd/"

# INITIAL CONDITIONS
const initial_cond = "rest"   # "rest" or "ncfile"
const initpath = "/network/aopp/chaos/pred/kloewer/julsdata/ssthr/"

const init_run_id = 2           # only for starting from ncfile

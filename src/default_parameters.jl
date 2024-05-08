@with_kw mutable struct Parameter

    T=Float32                 # number format

    Tprog=T                   # number format for prognostic variables
    Tcomm=Tprog               # number format for ghost-point copies
    Tini=Tprog                # number format to reduce precision for initial conditions

    # DOMAIN RESOLUTION AND RATIO
    nx::Int=100                         # number of grid cells in x-direction
    Lx::Float64=4000e3                  # length of the domain in x-direction [m]
    L_ratio::Float64=2                  # Domain aspect ratio of Lx/Ly

    # PHYSICAL CONSTANTS
    g::Float64 = 0.1                    # gravitational acceleration [m/s]
    H::Float64 = 500.                   # layer thickness at rest [m]
    ρ::Float64 = 1e3                    # water density [kg/m^3]
    ϕ::Float64 = 45.                    # central latitude of the domain (for coriolis) [°]
    ω::Float64 = 2π/(24*3600)           # Earth's angular frequency [s^-1]
    R::Float64 = 6.371e6                # Earth's radius [m]

    # SCALE
    scale::Float64=2^6                  # multiplicative scale for the momentum equations u,v
    scale_sst::Float64=2^15             # multiplicative scale for sst

    # WIND FORCING OPTIONS
    wind_forcing_x::String="shear"      # "channel", "double_gyre", "shear","constant" or "none"
    wind_forcing_y::String="constant"   # "channel", "double_gyre", "shear","constant" or "none"
    Fx0::Float64=0.12                   # wind stress strength [Pa] in x-direction
    Fy0::Float64=0.0                    # wind stress strength [Pa] in y-direction
    seasonal_wind_x::Bool=true          # Change the wind stress with a sine of frequency ωFx,ωFy
    seasonal_wind_y::Bool=false         # same for y-component
    ωFx::Float64=2                      # frequency [1/year] for x component
    ωFy::Float64=2                      # frequency [1/year] for y component

    # BOTTOM TOPOGRAPHY OPTIONS
    topography::String="ridges"         # "ridge", "seamount", "flat", "ridges", "bathtub"
    topo_ridges_positions::Vector{Float64} = [0.05,0.25,0.45,0.9]
    topo_height::Float64 = 100.         # height of seamount [m]
    topo_width::Float64 = 300e3         # horizontal scale [m] of the seamount

    # SURFACE RELAXATION
    surface_relax::Bool=false           # yes?
    t_relax::Float64 = 100.             # time scale of the relaxation [days]
    η_refh::Float64 = 5.                # height difference [m] of the interface relaxation profile
    η_refw::Float64 = 50e3              # width [m] of the tangent used for the interface relaxation

    # SURFACE FORCING
    surface_forcing::Bool=false         # yes?
    ωFη::Float64 = 1.0                  # frequency [1/year] for surfance forcing
    A::Float64 = 3e-5                   # Amplitude [m/s]
    ϕk::Float64 = ϕ                     # Central latitude of Kelvin wave pumping
    wk::Float64 = 10e3                  # width [m] in y of Gaussian used for surface forcing

    # TIME STEPPING OPTIONS
    time_scheme::String="RK"            # Runge-Kutta ("RK") or strong-stability preserving RK
                                        # "SSPRK2","SSPRK3","4SSPRK3"
    RKo::Int=4                          # Order of the RK time stepping scheme (2, 3 or 4)
    RKs::Int=3                          # Number of stages for SSPRK2
    RKn::Int=5                          # n^2 = s = Number of stages  for SSPRK3
    cfl::Float64 = 0.9                  # CFL number (1.0 recommended for RK4, 0.6 for RK3)
    Ndays::Float64 = 200.0              # number of days to integrate for
    nstep_diff::Int=1                   # diffusive part every nstep_diff time steps.
    nstep_advcor::Int=0                 # advection and coriolis update every nstep_advcor time steps.
                                        # 0 means it is included in every RK4 substep
    compensated::Bool=false             # Compensated summation in the time integration?

    # BOUNDARY CONDITION OPTIONS
    bc::String="periodic"               # "periodic" or anything else for nonperiodic
    α::Float64 = 2                      # lateral boundary condition parameter
                                        # 0 free-slip, 0<α<2 partial-slip, 2 no-slip

    # PARAMETERS FOR ADJOINT METHOD
    data_steps::StepRange{Int,Int} = 0:1:0      # Timesteps where data exists
    data::Array{Float32, 1} = [0.]              # model data
    J::Float64 = 0.                             # Placeholder for cost function evaluation
    j::Int = 1                                  # For keeping track of the entry in data

    # CHECKPOINTING VARIABLES
    i::Int = 0                                  # Placeholder for current timestep, needed for Checkpointing.jl

    # MOMENTUM ADVECTION OPTIONS
    adv_scheme::String="ArakawaHsu"     # "Sadourny" or "ArakawaHsu"
    dynamics::String="nonlinear"        # "linear" or "nonlinear"

    # BOTTOM FRICTION OPTIONS
    bottom_drag::String="none"          # "linear", "quadratic" or "none"
    cD::Float64 = 1e-5                  # bottom drag coefficient [dimensionless] for quadratic
    τD::Float64 = 300.                  # bottom drag coefficient [days] for linear

    # DIFFUSION OPTIONS
    diffusion::String="constant"        # "Smagorinsky" or "constant", biharmonic in both cases
    νB::Float64 = 500.0                 # [m^2/s] scaling constant for constant biharmonic diffusion
    cSmag::Float64 = 0.15               # Smagorinsky coefficient [dimensionless]

    # TRACER ADVECTION
    tracer_advection::Bool=true         # yes?
    tracer_relaxation::Bool=true        # yes?
    tracer_consumption::Bool=false      # yes?
    sst_initial::String="waves"         # "west", "south", "linear", "waves","rect", "flat" or "restart"
    sst_rect_coords::Array{Float64,1}=[0.,0.15,0.,1.0]
                                        # (x0,x1,y0,y1) are the size of the rectangle in [0,1]
    Uadv::Float64 = 0.2                 # Velocity scale [m/s] for tracer advection
    SSTmax::Float64 = 1.                # tracer (sea surface temperature) max for initial conditions
    SSTmin::Float64 = -1.               # tracer (sea surface temperature) min for initial conditions
    τSST::Float64 = 100                 # tracer restoring time scale [days]
    jSST::Float64 = 365                 # tracer consumption [days]
    SSTw::Float64 = 5e5                 # width [m] of the tangent used for the IC and interface relaxation
    SSTϕ::Float64 = 0.5                 # latitude/longitude fraction ∈ [0,1] of sst edge
    SSTwaves_ny::Float64 = 4            # wave crests/troughs in y
    SSTwaves_nx::Float64 = SSTwaves_ny*L_ratio  # wave crests/troughs in x
    SSTwaves_p::Float64 = 1/2           # power for rectangles (p<1)/smootheness(p>=1) of waves

    # OUTPUT OPTIONS
    output::Bool=false                  # netcdf output?
    output_vars::Array{String,1}=["u","v","η","sst"]  # which variables to output? "du","dv","dη","H","ζ" also allowed
    output_dt::Float64 = 24             # output time step [hours]
    outpath::String=pwd()               # path to output folder
    compression_level::Int=3            # compression level
    return_time::Bool=false             # return time of simulation of progn vars?

    # INITIAL CONDITIONS
    initial_cond::String="rest"         # "rest" or "ncfile" for restart from file
    initpath::String=outpath            # folder where to pick the restart files from
    init_run_id::Int=0                  # run id for restart from run number
    init_starti::Int=-1                 # timestep to start from (-1 meaning last)
    get_id_mode::String="continue"      # How to determine the run id: "continue" or "fill"
    run_id::Int=-1                      # Output with a specific run id
    init_interpolation::Bool=true       # Interpolate the initial conditions in case grids don't match?

    # ASSERT - CHECK THAT THE INPUT PARAMETERS MAKE SENSE
    @assert all((nx,Lx,L_ratio) .> 0.)  "nx, Lx, L_ratio have to be >0"
    @assert all((g,H,ρ,ω,R) .> 0.)      "g,H,ρ,ω,R have to be >0"
    @assert ϕ <= 90.0 && ϕ >= -90.0     "ϕ has to be in (-90,90), $ϕ given."
    @assert wind_forcing_x in ["channel","double_gyre","shear","constant","none"] "Wind forcing '$wind_forcing_x' unsupported"
    @assert wind_forcing_y in ["channel","double_gyre","shear","constant","none"] "Wind forcing '$wind_forcing_y' unsupported"
    @assert topography in ["ridge","seamount","flat","ridges"] "Topography '$topography' unsupported"
    @assert topo_width > 0.0    "topo_width has to be >0, $topo_width given."
    @assert t_relax > 0.0       "t_relax has to be >0, $t_relax given."
    @assert η_refw > 0.0        "η_refw has to be >0, $η_refw given."
    @assert time_scheme in ["RK","SSPRK2","SSPRK3","4SSPRK3"] "Time scheme $time_scheme unsupported."
    @assert RKo in [2,3,4]        "RKo has to be 2,3 or 4; $RKo given."
    @assert RKs > 1               "RKs has to be >= 2; $RKs given."
    @assert Ndays > 0.0         "Ndays has to be >0, $Ndays given."
    @assert nstep_diff > 0      "nstep_diff has to be >0, $nstep_diff given."
    @assert nstep_advcor >= 0   "nstep_advcor has to be >=0, $nstep_advcor given."
    @assert bc in ["periodic","nonperiodic"]    "boundary condition '$bc' unsupported."
    @assert α >= 0.0 && α <= 2.0    "Tangential boundary condition α has to be in (0,2), $α given."
    @assert adv_scheme in ["Sadourny","ArakawaHsu"] "Advection scheme '$adv_scheme' unsupported"
    @assert dynamics in ["linear","nonlinear"]  "Dynamics '$dynamics' unsupported."
    @assert bottom_drag in ["quadratic","linear","none"] "Bottom drag '$bottom_drag' unsupported."
    @assert cD >= 0.0    "Bottom drag coefficient cD has to be >=0, $cD given."
    @assert τD >= 0.0    "Bottom drag coefficient τD has to be >=0, $τD given."
    @assert diffusion in ["Smagorinsky", "constant"] "Diffusion '$diffusion' unsupported."
    @assert νB > 0.0     "Diffusion scaling constant νB has to be >0, $νB given."
    @assert cSmag > 0.0  "Smagorinsky coefficient cSmag has to be >0, $cSmag given."
    @assert Uadv > 0.0   "Advection velocity scale Uadv has to be >0, $Uadv given."
    @assert output_dt > 0   "Output time step has to be >0, $output_dt given."
    @assert initial_cond in ["rest", "ncfile"] "Initial conditions '$initial_cond' unsupported."
    @assert init_run_id >= 0 "Initial condition run id, init_run_id, has to be >= 0, $init_run_id given."
    @assert init_starti > 0 || init_starti == -1 "Start index, init_starti, has to be >0 || -1, $init_starti given."
    @assert get_id_mode in ["continue","fill","specific"] "get_id_mode $get_id_mode unsupported."
end

"""
Creates a Parameter struct with following options and default values
    T=Float32                 # number format

    Tprog=T                   # number format for prognostic variables
    Tcomm=Tprog               # number format for ghost-point copies
    Tini=Tprog                # number format to reduce precision for initial conditions

    # DOMAIN RESOLUTION AND RATIO
    nx::Int=100                         # number of grid cells in x-direction
    Lx::Float64=4000e3                  # length of the domain in x-direction [m]
    L_ratio::Float64=2                  # Domain aspect ratio of Lx/Ly

    # PHYSICAL CONSTANTS
    g::Float64 = 0.1                    # gravitational acceleration [m/s]
    H::Float64 = 500.                   # layer thickness at rest [m]
    ρ::Float64 = 1e3                    # water density [kg/m^3]
    ϕ::Float64 = 45.                    # central latitude of the domain (for coriolis) [°]
    ω::Float64 = 2π/(24*3600)           # Earth's angular frequency [s^-1]
    R::Float64 = 6.371e6                # Earth's radius [m]

    # SCALE
    scale::Float64=2^6                  # multiplicative scale for the momentum equations u,v
    scale_sst::Float64=2^15             # multiplicative scale for sst

    # WIND FORCING OPTIONS
    wind_forcing_x::String="shear"      # "channel", "double_gyre", "shear","constant" or "none"
    wind_forcing_y::String="constant"   # "channel", "double_gyre", "shear","constant" or "none"
    Fx0::Float64=0.12                   # wind stress strength [Pa] in x-direction
    Fy0::Float64=0.0                    # wind stress strength [Pa] in y-direction
    seasonal_wind_x::Bool=true          # Change the wind stress with a sine of frequency ωFx,ωFy
    seasonal_wind_y::Bool=false         # same for y-component
    ωFx::Float64=2                      # frequency [1/year] for x component
    ωFy::Float64=2                      # frequency [1/year] for y component

    # BOTTOM TOPOGRAPHY OPTIONS
    topography::String="ridges"         # "ridge", "seamount", "flat", "ridges", "bathtub"
    topo_ridges_positions::Vector{Float64} = [0.05,0.25,0.45,0.9]
    topo_height::Float64 = 100.         # height of seamount [m]
    topo_width::Float64 = 300e3         # horizontal scale [m] of the seamount

    # SURFACE RELAXATION
    surface_relax::Bool=false           # yes?
    t_relax::Float64 = 100.             # time scale of the relaxation [days]
    η_refh::Float64 = 5.                # height difference [m] of the interface relaxation profile
    η_refw::Float64 = 50e3              # width [m] of the tangent used for the interface relaxation

    # SURFACE FORCING
    surface_forcing::Bool=false         # yes?
    ωFη::Float64 = 1.0                  # frequency [1/year] for surfance forcing
    A::Float64 = 3e-5                   # Amplitude [m/s]
    ϕk::Float64 = ϕ                     # Central latitude of Kelvin wave pumping
    wk::Float64 = 10e3                  # width [m] in y of Gaussian used for surface forcing

    # TIME STEPPING OPTIONS
    time_scheme::String="RK"            # Runge-Kutta ("RK") or strong-stability preserving RK
                                        # "SSPRK2","SSPRK3","4SSPRK3"
    RKo::Int=4                          # Order of the RK time stepping scheme (2, 3 or 4)
    RKs::Int=3                          # Number of stages for SSPRK2
    RKn::Int=5                          # n^2 = s = Number of stages  for SSPRK3
    cfl::Float64 = 0.9                  # CFL number (1.0 recommended for RK4, 0.6 for RK3)
    Ndays::Float64 = 200.0              # number of days to integrate for
    nstep_diff::Int=1                   # diffusive part every nstep_diff time steps.
    nstep_advcor::Int=0                 # advection and coriolis update every nstep_advcor time steps.
                                        # 0 means it is included in every RK4 substep
    compensated::Bool=false             # Compensated summation in the time integration?

    # BOUNDARY CONDITION OPTIONS
    bc::String="periodic"               # "periodic" or anything else for nonperiodic
    α::Float64 = 2                      # lateral boundary condition parameter
                                        # 0 free-slip, 0<α<2 partial-slip, 2 no-slip

    # PARAMETERS FOR ADJOINT METHOD
    data_steps::StepRange{Int,Int} = 0:1:0      # Timesteps where data exists
    data::Array{Float32, 1} = [0.]              # model data
    J::Float64 = 0.                             # Placeholder for cost function evaluation
    j::Int = 1                                  # For keeping track of the entry in data

    # CHECKPOINTING VARIABLES
    i::Int = 0                                  # Placeholder for current timestep, needed for Checkpointing.jl

    # MOMENTUM ADVECTION OPTIONS
    adv_scheme::String="ArakawaHsu"     # "Sadourny" or "ArakawaHsu"
    dynamics::String="nonlinear"        # "linear" or "nonlinear"

    # BOTTOM FRICTION OPTIONS
    bottom_drag::String="none"          # "linear", "quadratic" or "none"
    cD::Float64 = 1e-5                  # bottom drag coefficient [dimensionless] for quadratic
    τD::Float64 = 300.                  # bottom drag coefficient [days] for linear

    # DIFFUSION OPTIONS
    diffusion::String="constant"        # "Smagorinsky" or "constant", biharmonic in both cases
    νB::Float64 = 500.0                 # [m^2/s] scaling constant for constant biharmonic diffusion
    cSmag::Float64 = 0.15               # Smagorinsky coefficient [dimensionless]

    # TRACER ADVECTION
    tracer_advection::Bool=true         # yes?
    tracer_relaxation::Bool=true        # yes?
    tracer_consumption::Bool=false      # yes?
    sst_initial::String="waves"         # "west", "south", "linear", "waves","rect", "flat" or "restart"
    sst_rect_coords::Array{Float64,1}=[0.,0.15,0.,1.0]
                                        # (x0,x1,y0,y1) are the size of the rectangle in [0,1]
    Uadv::Float64 = 0.2                 # Velocity scale [m/s] for tracer advection
    SSTmax::Float64 = 1.                # tracer (sea surface temperature) max for initial conditions
    SSTmin::Float64 = -1.               # tracer (sea surface temperature) min for initial conditions
    τSST::Float64 = 100                 # tracer restoring time scale [days]
    jSST::Float64 = 365                 # tracer consumption [days]
    SSTw::Float64 = 5e5                 # width [m] of the tangent used for the IC and interface relaxation
    SSTϕ::Float64 = 0.5                 # latitude/longitude fraction ∈ [0,1] of sst edge
    SSTwaves_ny::Float64 = 4            # wave crests/troughs in y
    SSTwaves_nx::Float64 = SSTwaves_ny*L_ratio  # wave crests/troughs in x
    SSTwaves_p::Float64 = 1/2           # power for rectangles (p<1)/smootheness(p>=1) of waves

    # OUTPUT OPTIONS
    output::Bool=false                  # netcdf output?
    output_vars::Array{String,1}=["u","v","η","sst"]  # which variables to output? "du","dv","dη","H","ζ" also allowed
    output_dt::Float64 = 24             # output time step [hours]
    outpath::String=pwd()               # path to output folder
    compression_level::Int=3            # compression level
    return_time::Bool=false             # return time of simulation of progn vars?

    # INITIAL CONDITIONS
    initial_cond::String="rest"         # "rest" or "ncfile" for restart from file
    initpath::String=outpath            # folder where to pick the restart files from
    init_run_id::Int=0                  # run id for restart from run number
    init_starti::Int=-1                 # timestep to start from (-1 meaning last)
    get_id_mode::String="continue"      # How to determine the run id: "continue" or "fill"
    run_id::Int=-1                      # Output with a specific run id
    init_interpolation::Bool=true       # Interpolate the initial conditions in case grids don't match?
"""
Parameter
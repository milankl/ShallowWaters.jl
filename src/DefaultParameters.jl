@with_kw struct Parameter

    T::DataType=Float32                 # number format

    # DOMAIN RESOLUTION AND RATIO
    nx::Int=100                         # number of grid cells in x-direction
    Lx::Real=2000e3                     # length of the domain in x-direction [m]
    L_ratio::Real=2                     # Domain aspect ratio of Lx/Ly

    # PHYSICAL CONSTANTS
    g::Real=10.                         # gravitational acceleration [m/s]
    H::Real=500.                        # layer thickness at rest [m]
    ρ::Real=1e3                         # water density [kg/m^3]
    ϕ::Real=15.                         # central latitue of the domain (for coriolis) [°]
    ω::Real=2π/(24*3600)                # Earth's angular frequency [s^-1]
    R::Real=6.371e6                     # Earth's radius [m]

    # WIND FORCING OPTIONS
    wind_forcing_x::String="channel"    # "channel", "double_gyre", "shear","constant" or "none"
    wind_forcing_y::String="constant"   # "channel", "double_gyre", "shear","constant" or "none"
    Fx0::Real=0.12                      # wind stress strength [Pa] in x-direction
    Fy0::Real=0.0                       # wind stress strength [Pa] in y-direction

    # BOTTOM TOPOGRAPHY OPTIONS
    topography::String="ridge"          # "ridge", "seamount", "flat", "ridges", "bathtub"
    topo_height::Real=10.               # height of seamount [m]
    topo_width::Real=300e3              # horizontal scale [m] of the seamount

    # SURFACE RELAXATION
    surface_relax::Bool=false           # yes?
    t_relax::Real=100.                  # time scale of the relaxation [days]
    η_refh::Real=5.                     # height difference [m] of the interface relaxation profile
    η_refw::Real=100e3                  # width [m] of the tangent used for the interface relaxation

    # SURFACE FORCING
    surface_forcing::Bool=false         # yes?
    ωyr::Real=1.0                       # (annual) frequency [1/year]
    A::Real=3e-5                        # Amplitude [m/s]

    # TIME STEPPING OPTIONS
    RKo::Int=4                          # Order of the RK time stepping scheme (3 or 4)
    cfl::Real=1.0                       # CFL number (1.0 recommended for RK4, 0.6 for RK3)
    Ndays::Real=10.0                    # number of days to integrate for
    nstep_diff::Int=1                   # diffusive part every nstep_diff time steps.
    nstep_advcor::Int=0                 # advection and coriolis update every nstep_advcor time steps.
                                        # 0 means it is included in every RK4 substep

    # BOUNDARY CONDITION OPTIONS
    bc::String="periodic"               # "periodic" or anything else for nonperiodic
    α::Real=2.                          # lateral boundary condition parameter
                                        # 0 free-slip, 0<α<2 partial-slip, 2 no-slip

    # MOMENTUM ADVECTION OPTIONS
    adv_scheme::String="ArakawaHsu"     # "Sadourny" or "ArakawaHsu"
    dynamics::String="nonlinear"        # "linear" or "nonlinear"

    # BOTTOM FRICTION OPTIONS
    bottom_drag::String="quadratic"     # "linear", "quadratic" or "none"
    cD::Real=1e-5                       # bottom drag coefficient [dimensionless] for quadratic
    τD::Real=300.                       # bottom drag coefficient [days] for linear

    # DIFFUSION OPTIONS
    diffusion::String="Smagorinsky"     # "Smagorinsky" or "Constant", biharmonic in both cases
    νB::Real=500.0                      # [m^2/s] scaling constant for Constant biharmonic diffusion
    cSmag::Real=0.15                    # Smagorinsky coefficient [dimensionless]

    # TRACER ADVECTION
    tracer_advection::Bool=false        # yes?
    tracer_relaxation::Bool=false       # yes?
    tracer_consumption::Bool=false      # yes?
    tracer_pumping::Bool=false          # yes?
    injection_region::String="west"     # "west", "south", "rect" or flat
    sst_initial::String="south"         # same here
    sstrestart::Bool=true               # start from previous sst file
    Uadv::Real=0.15                     # Velocity scale [m/s] for tracer advection
    SSTmax::Real=1.                     # tracer (sea surface temperature) max for restoring
    SSTmin::Real=0.                     # tracer (sea surface temperature) min for restoring
    τSST::Real=500.                     # tracer restoring time scale [days]
    jSST::Real=50*365.                  # tracer consumption [days]
    SST_λ0::Real=222e3                  # [m] transition position of relaxation timescale
    SST_λs::Real=111e3                  # [m] transition width of relaxation timescale
    SST_γ0::Real=8.35                   # [days] injection time scale
    SSTw::Real=1000e3                   # width [m] of the tangent used for the IC and interface relaxation
    SSTϕ::Real=0.5                      # latitude/longitude fraction ∈ [0,1] of sst edge

    # OUTPUT OPTIONS
    output::Bool=false                  # for nc output
    output_tend::Bool=false             # ouput for tendencies as well?
    output_diagn::Bool=false            # output diagnostic variables as well?

                                        # which prognostic variables to output?
    output_progn_vars::Array{String,1}=["u","v","eta","sst"]
                                        # which tendencies to output?
    output_tend_vars::Array{String,1}=["du","dv","deta","Bu","Bv","LLu1","LLu2","LLv1","LLv2",
                                "qhv","qhu","dpdx","dpdy","dUdx","dVdy"]
                                        # which diagnostic variables to output?
    output_diagn_vars::Array{String,1}=["q","p","dudx","dvdy","dudy","dvdx","Lu",
                                "Lv","xd","yd"]

    output_dt::Real=6                   # output time step [hours]
    outpath::String="data/"             # path to output folder

    # INITIAL CONDITIONS
    initial_cond::String="rest"         # "rest" or "ncfile" for restart from file
    initpath::String="data/"            # folder where to pick the restart files from
    init_run_id::Int=0                  # run id for restart from run number

    # ASSERT - CHECK THAT THE INPUT PARAMETERS MAKE SENSE
    @assert all((nx,Lx,L_ratio) .> 0.)
    @assert all((g,H,ρ,ω,R) .> 0.)
    @assert ϕ <= 90.0 && ϕ >= -90.0
    #TODO more of that

end

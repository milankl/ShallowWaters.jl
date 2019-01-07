"""Start from rest or load initial conditions from netCDF file."""
function initial_conditions(starti=-1)
    # initialise the state matrices u,v,η and set their initial conditions

    if initial_cond == "rest"
        # sub domain sub_dom is the entire domain for a single process (no domain decomposition)
        # u = zeros(Numtype,sub_dom["nux"],sub_dom["nuy"])
        # v = zeros(Numtype,sub_dom["nvx"],sub_dom["nvy"])
        # η = zeros(Numtype,sub_dom["nx"],sub_dom["ny"])
        u = zeros(Numtype,nux,nuy)
        v = zeros(Numtype,nvx,nvy)
        η = zeros(Numtype,nx,ny)
        sst = sst_initial()

    elseif initial_cond == "ncfile"

        #TODO for domain decomposition: slice ncfile for different processors
        inirunpath = initpath*"run"*@sprintf("%04d",init_run_id)*"/"

        # take starti time step from existing netcdf files
        ncu = NetCDF.open(inirunpath*"u.nc")

        if starti == -1
            # replace -1 with length of time dimension
            starti = size(ncu.vars["t"])[end]
        end

        u = ncu.vars["u"][:,:,starti]
        NetCDF.close(ncu)

        ncv = NetCDF.open(inirunpath*"v.nc")
        v = ncv.vars["v"][:,:,starti]
        NetCDF.close(ncv)

        ncη = NetCDF.open(inirunpath*"eta.nc")
        η = ncη.vars["eta"][:,:,starti]
        NetCDF.close(ncη)

        if sstrestart
            ncsst = NetCDF.open(inirunpath*"sst.nc")
            sst = ncsst.vars["sst"][:,:,starti]
            NetCDF.close(ncsst)

            sst = Numtpe.(reshape(sst,size(sst)[1:2]))
        else
            sst = sst_initial()
        end

        # remove singleton time dimension
        # and convert from Float32 to Numtype
        u = Numtype.(reshape(u,size(u)[1:2]))
        v = Numtype.(reshape(v,size(v)[1:2]))
        η = Numtype.(reshape(η,size(η)[1:2]))
    end

    return u,v,η,sst
end

"""Zonally constant hyperbolic tangent initial conditions for the tracer determined by SSTmax, SSTmin, SSTw, SSTϕ."""
function sst_south()
    xx_T,yy_T = meshgrid(x_T,y_T)
    sst = (SSTmin+SSTmax)/2 .+ tanh.(2π*(Ly/(4*SSTw))*(yy_T/Ly .- SSTϕ))*(SSTmin-SSTmax)/2
    return Numtype.(sst)
end

"""Meriodionally constant hyperbolic tangent initial conditions for the tracer determined by SSTmax, SSTmin, SSTw, SSTϕ"""
function sst_west()
    xx_T,yy_T = meshgrid(x_T,y_T)
    sst = (SSTmin+SSTmax)/2 .+ tanh.(2π*(Lx/(4*SSTw))*(xx_T/Lx .- SSTϕ))*(SSTmin-SSTmax)/2
    return Numtype.(sst)
end

"""Initial conditions as a rectangle with SSTmax inside, SSTmin outside."""
function sst_rect()

    x0,x1 = 0.05,0.45     # left right border in [0,1]
    y0,y1 = 0.,1.0     # north south border in [0,1]

    xx_T,yy_T = meshgrid(x_T,y_T)

    sst = fill(SSTmin,size(xx_T))
    inside = (xx_T/Lx .> x0) .* (xx_T/Lx .< x1) .* (yy_T/Ly .> y0) .* (yy_T/Ly .< y1)
    sst[inside] .= SSTmax
    return Numtype.(sst)
end

if injection_region == "south"
    sst_initial = sst_south
elseif injection_region == "west"
    sst_initial = sst_west
elseif injection_region == "rect"
    sst_initial = sst_rect
end

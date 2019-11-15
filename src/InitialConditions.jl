"""
    P = ProgVars{T}(u,v,η,sst)

Struct containing the prognostic variables u,v,η and sst.
"""
struct PrognosticVars{T<:AbstractFloat}
    u::Array{T,2}           # u-velocity
    v::Array{T,2}           # v-velocity
    η::Array{T,2}           # sea surface height / interface displacement
    sst::Array{T,2}         # tracer / sea surface temperature
end

function initial_conditions(::Type{T},S::ModelSetup) where {T<:AbstractFloat}

    ## PROGNOSTIC VARIABLES U,V,η

    @unpack nux,nuy,nvx,nvy,nx,ny = S.grid
    @unpack initial_cond = S.parameters

    if intitial_cond == "rest"

        u = zeros(T,nux,nuy)
        v = zeros(T,nvx,nvy)
        η = zeros(T,nx,ny)

    elseif initial_cond == "ncfile"

        @unpack initpath,init_run_id,init_starti = S.parameters

        inirunpath = initpath*"run"*@sprintf("%04d",init_run_id)*"/"

        # take starti time step from existing netcdf files
        ncu = NetCDF.open(inirunpath*"u.nc")

        if init_starti == -1    # replace -1 with length of time dimension
            init_starti = size(ncu.vars["t"])[end]
        end

        u = ncu.vars["u"][:,:,init_starti]
        NetCDF.close(ncu)

        ncv = NetCDF.open(inirunpath*"v.nc")
        v = ncv.vars["v"][:,:,init_starti]
        NetCDF.close(ncv)

        ncη = NetCDF.open(inirunpath*"eta.nc")
        η = ncη.vars["eta"][:,:,init_starti]
        NetCDF.close(ncη)

        # remove singleton time dimension
        # and convert from Float32 to Numtype
        u = T.(reshape(u,size(u)[1:2]))
        v = T.(reshape(v,size(v)[1:2]))
        η = T.(reshape(η,size(η)[1:2]))

    end

    ## SST

    @unpack SSTmin, SSTmax, SSTw, SSTϕ = S.parameters
    @unpack sst_initial = S.parameters
    @unpack x_T, y_T, Ly = S.grid

    xx_T,yy_T = meshgrid(x_T,y_T)

    if sst_initial == "south"
        sst = (SSTmin+SSTmax)/2 .+ tanh.(2π*(Ly/(4*SSTw))*(yy_T/Ly .- SSTϕ))*(SSTmin-SSTmax)/2
    elseif sst_initial == "west"
        sst = (SSTmin+SSTmax)/2 .+ tanh.(2π*(Lx/(4*SSTw))*(xx_T/Lx .- SSTϕ))*(SSTmin-SSTmax)/2
    elseif sst_initial == "flat"
        sst = fill(SSTmin,size(xx_T))
    elseif sst_initial == "rect"
        @unpack sst_rect_coords = S.parameters
        x0,x1,y0,y1 = sst_rect_coords

        sst = fill(SSTmin,size(xx_T))
        inside = (xx_T/Lx .> x0) .* (xx_T/Lx .< x1) .* (yy_T/Ly .> y0) .* (yy_T/Ly .< y1)
        sst[inside] .= SSTmax
    end

    if initial_cond == "ncfile" && sst_initial == "restart"
        ncsst = NetCDF.open(inirunpath*"sst.nc")
        sst = ncsst.vars["sst"][:,:,init_starti]
        NetCDF.close(ncsst)

        sst = reshape(sst,size(sst)[1:2])
    end

    # Convert to number format T
    sst = T.(sst)

    u,v,η,sst = add_halo(u,v,η,sst,S)

    return PrognosticVars{T}(u,v,η,sst)
end

# if injection_region == "south"
#     sst_inj_region = sst_south
# elseif injection_region == "west"
#     sst_inj_region = sst_west
# elseif injection_region == "rect"
#     sst_inj_region = sst_rect
# end

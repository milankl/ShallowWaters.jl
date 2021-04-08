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

"""Zero generator function for Grid G as argument."""
function PrognosticVars{T}(G::Grid) where {T<:AbstractFloat}
    @unpack nux,nuy,nvx,nvy,nx,ny = G
    @unpack halo,haloη,halosstx,halossty = G

    u = zeros(T,nux+2halo,nuy+2halo)
    v = zeros(T,nvx+2halo,nvy+2halo)
    η = zeros(T,nx+2haloη,ny+2haloη)
    sst = zeros(T,nx+2halosstx,ny+2halossty)

    return PrognosticVars{T}(u,v,η,sst)
end

function initial_conditions(::Type{T},S::ModelSetup) where {T<:AbstractFloat}

    ## PROGNOSTIC VARIABLES U,V,η
    @unpack nux,nuy,nvx,nvy,nx,ny = S.grid
    @unpack initial_cond = S.parameters
    @unpack Tini = S.parameters

    if initial_cond == "rest"

        u = zeros(T,nux,nuy)
        v = zeros(T,nvx,nvy)
        η = zeros(T,nx,ny)

    elseif initial_cond == "ncfile"

        @unpack initpath,init_run_id,init_starti = S.parameters
        @unpack init_interpolation = S.parameters
        @unpack nx,ny = S.grid

        inirunpath = joinpath(initpath,"run"*@sprintf("%04d",init_run_id))

        # take starti time step from existing netcdf files
        ncstring = joinpath(inirunpath,"u.nc")
        ncu = NetCDF.open(ncstring)

        if init_starti == -1    # replace -1 with length of time dimension
            init_starti = size(ncu.vars["t"])[end]
        end

        u = ncu.vars["u"][:,:,init_starti]

        ncv = NetCDF.open(joinpath(inirunpath,"v.nc"))
        v = ncv.vars["v"][:,:,init_starti]

        ncη = NetCDF.open(joinpath(inirunpath,"eta.nc"))
        η = ncη.vars["eta"][:,:,init_starti]

        # remove singleton time dimension
        u = reshape(u,size(u)[1:2])
        v = reshape(v,size(v)[1:2])
        η = reshape(η,size(η)[1:2])

        nx_old,ny_old = size(η)

        if (nx_old,ny_old) != (nx,ny)
            if init_interpolation           # interpolation in case the grids don't match

                # old grids
                x_T = collect(0.5:nx_old-0.5)
                y_T = collect(0.5:ny_old-0.5)

                # assuming periodic BC for now #TODO make flexible
                x_u = collect(0:nx_old-1)
                y_u = y_T

                x_v = x_T
                y_v = collect(1:ny_old-1)

                # set up interpolation functions
                u_itp = interpolate((x_u,y_u),u,Gridded(Linear()))
                v_itp = interpolate((x_v,y_v),v,Gridded(Linear()))
                η_itp = interpolate((x_T,y_T),η,Gridded(Linear()))

                #TODO in case of interpolation on larger grids
                #TODO make BC adapting to actual BCs used.
                u_etp = extrapolate(u_itp,(Flat(),Flat()))
                v_etp = extrapolate(v_itp,(Flat(),Flat()))
                η_etp = extrapolate(η_itp,(Flat(),Flat()))

                # new grids
                Δx = nx_old/nx
                Δy = ny_old/ny

                x_T_new = collect(Δx/2:Δx:nx_old-Δx/2)
                y_T_new = collect(Δy/2:Δy:ny_old-Δy/2)

                x_u_new = collect(0:Δx:nx_old-Δx)
                y_u_new = y_T_new

                x_v_new = x_T_new
                y_v_new = collect(Δy:Δy:ny_old-Δy)

                # retrieve values and overwrite existing arrays
                u = u_etp(x_u_new,y_u_new)
                v = v_etp(x_v_new,y_v_new)
                η = η_etp(x_T_new,y_T_new)

            else
                throw(error("Grid size $((nx,ny)) doesn't match initial conditions on a $(size(η)) grid."))
            end
        end
    end

    ## SST

    @unpack SSTmin, SSTmax, SSTw, SSTϕ = S.parameters
    @unpack sst_initial,scale = S.parameters
    @unpack x_T,y_T,Lx,Ly = S.grid

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
        ncsst = NetCDF.open(joinpath(inirunpath,"sst.nc"))
        sst = ncsst.vars["sst"][:,:,init_starti]

        sst = reshape(sst,size(sst)[1:2])
    end

    # Convert to number format T
    # allow for comparable initial conditions via Tini
    sst = T.(Tini.(sst))
    u = T.(Tini.(u))
    v = T.(Tini.(v))
    η = T.(Tini.(η))

    #TODO SST INTERPOLATION
    u,v,η,sst = add_halo(u,v,η,sst,S)

    return PrognosticVars{T}(u,v,η,sst)
end

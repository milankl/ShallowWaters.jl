struct NcFiles
    u::Union{NcFile,Nothing}        # zonal velocity
    v::Union{NcFile,Nothing}        # meridional velocity
    η::Union{NcFile,Nothing}        # sea surface height anomaly
    sst::Union{NcFile,Nothing}      # tracer / sea surface temperature
    q::Union{NcFile,Nothing}        # potentival vorticity
    ζ::Union{NcFile,Nothing}        # relative vorticity
    iout::Array{Integer,1}          # output index, stored in array for mutability
end

"""Generator function for "empty" NcFiles struct."""
NcFiles(x::Nothing) = NcFiles(x,x,x,x,x,x,[0])

"""Generator function for NcFiles struct, creating the underlying netCDF files."""
function NcFiles(feedback::Feedback,S::ModelSetup)

    if S.parameters.output

        @unpack output_vars = S.parameters
        @unpack x_u,x_v,x_T,x_q = S.grid
        @unpack y_u,y_v,y_T,y_q = S.grid

        @unpack run_id,runpath = feedback

        ncu = if "u" in output_vars nc_create(x_u,y_u,"u",runpath,"m/s","zonal velocity") else nothing end
        ncv = if "v" in output_vars nc_create(x_v,y_v,"v",runpath,"m/s","meridional velocity") else nothing end
        ncη = if "η" in output_vars nc_create(x_T,y_T,"eta",runpath,"m","sea surface height") else nothing end
        ncsst = if "sst" in output_vars nc_create(x_T,y_T,"sst",runpath,"1","sea surface temperature") else nothing end
        ncq = if "q" in output_vars nc_create(x_q,y_q,"q",runpath,"1/(ms)","potential vorticity") else nothing end
        ncζ = if "ζ" in output_vars nc_create(x_q,y_q,"relvort",runpath,"1","relative vorticity") else nothing end

        for nc in (ncu,ncv,ncη,ncsst,ncq,ncζ)
            if nc != nothing
                NetCDF.putatt(nc,"t",Dict("units"=>"s","long_name"=>"time"))
                NetCDF.putatt(nc,"x",Dict("units"=>"m","long_name"=>"zonal coordinate"))
                NetCDF.putatt(nc,"y",Dict("units"=>"m","long_name"=>"meridional coordinate"))
            end
        end

        return NcFiles(ncu,ncv,ncη,ncsst,ncq,ncζ,[0])
    else
        return NcFiles(nothing)
    end
end

"""Creates a netCDF file based on grid vectors x,y the variable name, its path, unit and long_name."""
function nc_create( x::Array{T,1},
                    y::Array{T,1},
                    name::String,
                    path::String,
                    unit::String,
                    long_name::String) where {T<:Real}

    xdim = NcDim("x",length(x),values=x)
    ydim = NcDim("y",length(y),values=y)
    tdim = NcDim("t",0,unlimited=true)

    var = NcVar(name,[xdim,ydim,tdim],t=Float32)
    tvar = NcVar("t",tdim,t=Int32)

    nc = NetCDF.create(joinpath(path,name*".nc"),[var,tvar],mode=NC_NETCDF4)
    NetCDF.putatt(nc,name,Dict("units"=>unit,"long_name"=>long_name))
    return nc
end

function output_nc!(i::Int,
                    ncs::NcFiles,
                    Prog::PrognosticVars,
                    Diag::DiagnosticVars,
                    S::ModelSetup)

    @unpack nout = S.grid
    @unpack output = S.parameters

    if i % nout == 0 && output

        iout = ncs.iout[1]+1        # unpack and increase
        ncs.iout[1] = iout          # pack

        @unpack halo,haloη,halosstx,halossty = S.grid
        @unpack q,dvdx,dudy = Diag.Vorticity
        @unpack f_q,ep,dtint = S.grid

        # CUT OFF HALOS
        @views u = Prog.u[halo+1:end-halo,halo+1:end-halo]
        @views v = Prog.v[halo+1:end-halo,halo+1:end-halo]
        @views η = Prog.η[haloη+1:end-haloη,haloη+1:end-haloη]
        @views sst = Prog.sst[halosstx+1:end-halosstx,halossty+1:end-halossty]
        @views ζ = (dvdx[2:end-1,2:end-1]-dudy[2+ep:end-1,2:end-1])./abs.(f_q)

        # WRITING THE VARIABLES
        if ncs.u != nothing
            NetCDF.putvar(ncs.u,"u",Float32.(u),start=[1,1,iout],count=[-1,-1,1])
        end
        if ncs.v != nothing
            NetCDF.putvar(ncs.v,"v",Float32.(v),start=[1,1,iout],count=[-1,-1,1])
        end
        if ncs.η != nothing
            NetCDF.putvar(ncs.η,"eta",Float32.(η),start=[1,1,iout],count=[-1,-1,1])
        end
        if ncs.sst != nothing
            NetCDF.putvar(ncs.sst,"sst",Float32.(sst),start=[1,1,iout],count=[-1,-1,1])
        end
        if ncs.q != nothing
            NetCDF.putvar(ncs.q,"q",Float32.(q),start=[1,1,iout],count=[-1,-1,1])
        end
        if ncs.ζ != nothing
            NetCDF.putvar(ncs.ζ,"relvort",Float32.(ζ),start=[1,1,iout],count=[-1,-1,1])
        end

        # WRITING THE TIME
        for nc in (ncs.u,ncs.v,ncs.η,ncs.sst,ncs.q,ncs.ζ)
            if nc != nothing
                NetCDF.putvar(nc,"t",Int64[i*dtint],start=[iout])
                NetCDF.sync(nc)         # sync to view netcdf while model is still running
            end
        end
    end
end

"""Closes netCDF and progress.txt files."""
function output_close!(ncs::NcFiles,feedback::Feedback,S::ModelSetup)

    @unpack output = S.parameters

    if output
        for nc in (ncs.u,ncs.v,ncs.η,ncs.sst,ncs.q,ncs.ζ)
            if nc != nothing
                NetCDF.close(nc)
            end
        end
        println("All data stored.")
        write(feedback.progress_txt,"All data stored.")
        close(feedback.progress_txt)
    end
end

"""Finds the first gap in a list of integers."""
function gap(a::Array{Int,1})
    try
        return minimum([i for i in minimum(a):maximum(a) if ~(i in a)])
    catch
        return maximum(a)+1
    end
end

"""Checks output folders to determine a 4-digit run id number."""
function get_run_id_path(   S::ModelSetup,
                            order::String="continue",
                            run_id::Union{Int,Nothing}=nothing)

    @unpack output,outpath = S.parameters

    if output
        runlist = filter(x->startswith(x,"run"),readdir(outpath))
        existing_runs = [parse(Int,id[4:end]) for id in runlist]
        if length(existing_runs) == 0           # if no runfolder exists yet
            runpath = outpath*"run0000/"
            mkdir(runpath)
            return 0,runpath
        else                                    # create next folder
            if order == "fill"  # find the smallest gap in runfolders
                run_id = gap(existing_runs)
                runpath = outpath*"run"*@sprintf("%04d",run_id)*"/"
                mkdir(runpath)

            elseif order == "specific" # specify the run_id as input argument
                runpath = outpath*"run"*@sprintf("%04d",run_id)*"/"
                try # create folder if not existent
                    mkdir(runpath)
                catch # else rm folder and create new one
                    rm(runpath,recursive=true)
                    mkdir(runpath)
                end

            elseif order == "continue" # find largest folder and count one up
                run_id = maximum(existing_runs)+1
                runpath = outpath*"run"*@sprintf("%04d",run_id)*"/"
                mkdir(runpath)
            else
                throw(error("Order $order is not valid for get_run_id_path(), chose continue, specific or fill."))
            end
            return run_id,runpath
        end
    else
        return 0,"no runpath"
    end
end

# #TODO in ensemble mode, the .jl files might have changed since the start and do not correspond to what
# #TODO is actually executed!
# """Archives all .jl files of juls in the output folder to make runs reproducible."""
# function scripts_output()
#     if output #&& prank == 0
#         # copy all files in juls main folder
#         mkdir(runpath*"scripts")
#         for juliafile in filter(x->endswith(x,".jl"),readdir())
#             cp(juliafile,runpath*"scripts/"*juliafile)
#         end
#
#         # and also in the src folder
#         mkdir(runpath*"scripts/src")
#         for juliafile in filter(x->endswith(x,".jl"),readdir("src"))
#             cp("src/"*juliafile,runpath*"scripts/src/"*juliafile)
#         end
#     end
# end
#
# """Creates a dictionary with many parameter constants to be included in the nc files."""
# function output_dict()
#     # Attributes for nc
#     Dictu = Dict{String,Any}("description"=>"Data from shallow-water model juls.")
#     Dictu["details"] = "Cartesian coordinates, f or beta-plane, Arakawa C-grid"
#     Dictu["reference"] = "github.com/milankl/juls"
#
#     Dictu["nx"] = nx
#     Dictu["Lx"] = Lx
#     Dictu["L_ratio"] = L_ratio
#     Dictu["delta"] = Δ
#
#     Dictu["halo"] = halo
#     Dictu["haloeta"] = haloη
#     Dictu["halosstx"] = halosstx
#     Dictu["halossty"] = halossty
#
#     Dictu["g"] = gravity
#     Dictu["water_depth"] = water_depth
#     Dictu["phi"] = ϕ
#     Dictu["density"] = ρ
#
#     Dictu["wind_forcing_x"] = wind_forcing_x
#     Dictu["wind_forcing_y"] = wind_forcing_y
#     Dictu["Fx0"] = Fx0
#     Dictu["Fy0"] = Fy0
#
#     Dictu["topography_feature"] = topography_feature
#     Dictu["topofeat_height"] = topofeat_height
#     Dictu["topofeat_width"] = topofeat_width
#
#     Dictu["surface_forcing"] = string(surface_forcing)
#     Dictu["t_relax"] = t_relax
#     Dictu["eta_refh"] = η_refh
#     Dictu["η_refw"] = η_refw
#
#     Dictu["Numtype"] = string(Numtype)
#     Dictu["output_dt"] = output_dt
#     Dictu["nout"] = nout
#     Dictu["nadvstep"] = nadvstep
#     Dictu["nstep_diff"] = nstep_diff
#     Dictu["nstep_advcor"] = nstep_advcor
#
#     Dictu["RKo"] = RKo
#     Dictu["cfl"] = cfl
#     Dictu["Ndays"] = Ndays
#
#     Dictu["bc_x"] = bc_x
#     Dictu["lbc"] = lbc
#
#     Dictu["adv_scheme"] = adv_scheme
#     Dictu["dynamics"] = dynamics
#
#     Dictu["bottom_friction"] = bottom_friction
#     Dictu["drag"] = drag
#     Dictu["taudrag"] = τdrag
#
#     Dictu["diffusion"] = diffusion
#     Dictu["nuConst"] = ν_const
#     Dictu["c_smag"] = c_smag
#
#     Dictu["tracer_advcetion"] = string(tracer_advection)
#     Dictu["tracer_relaxation"] = string(tracer_relaxation)
#     Dictu["injection_region"] = injection_region
#     Dictu["sstrestart"] = string(sstrestart)
#     Dictu["Uadv"] = Uadv
#     Dictu["SSTmax"] = SSTmax
#     Dictu["SSTmin"] = SSTmin
#     Dictu["tauSST"] = τSST
#     Dictu["SSTw"] = SSTw
#     Dictu["SSTphi"] = SSTϕ
#
#     Dictu["initial_cond"] = initial_cond
#     Dictu["init_run_id"] = init_run_id
#     Dictu["initpath"] = initpath
#
#     return Dictu
# end

"""Initialises netCDF files for data output of the prognostic variables."""
function output_ini(u,v,η,sst,Bu,Bv,LLu1,LLu2,LLv1,LLv2)
    # only process with rank 0 defines the netCDF file
    if output #&& prank == 0

        all_output_progn_vars = ["u","v","eta","sst"]
        units_progn = ["m/s","m/s","m","degC"]
        longnames_progn = ["zonal velocity","meridional velocity","sea surface height","sea surface temperature"]

        all_output_diagn_vars = ["Bu","Bv","LLu1","LLu2","LLv1","LLv2"]
        units_diagn = ["m^2/s^2","m^2/s^2","m^2/s^2","m^2/s^2","m^2/s^2","m^2/s^2"]
        longnames_diagn = ["Bottom friction u-comp.","Bottom friction v-comp.",
                            "Diffusion u-comp. 1","Diffusion u-comp. 2",
                            "Diffusion v-comp. 1","Diffusion v-comp. 2"]

        allx = (x_u,x_v,x_T,x_q)    # collect all grids
        ally = (y_u,y_v,y_T,y_q)
        grids_progn = [1,2,3,3]           # for easy access per index
        grids_diagn = [1,2,1,1,2,2]

        ncs_progn = Array{Any,1}(zeros(Int,length(all_output_progn_vars)))
        ncs_diagn = Array{Any,1}(zeros(Int,length(all_output_diagn_vars)))

        # loop over all outputtable variables
        # PROGNOSTIC VARIABLES
        for (ivarout,outvar) in enumerate(all_output_progn_vars)
            if outvar in output_progn_vars    # check whether output is desired (specified in parameters.jl)
                ncs_progn[ivarout] = nccreate(allx[grids_progn[ivarout]],ally[grids_progn[ivarout]],
                                outvar,runpath,units_progn[ivarout],longnames_progn[ivarout])
            end
        end

        # DIAGNOSTIC VARIABLES
        if output_diagn
            for (ivarout,outvar) in enumerate(all_output_diagn_vars)
                if outvar in output_diagn_vars    # check whether output is desired (specified in parameters.jl)
                    ncs_diagn[ivarout] = nccreate(allx[grids_diagn[ivarout]],ally[grids_diagn[ivarout]],
                                    outvar,runpath,units_diagn[ivarout],longnames_diagn[ivarout])
                end
            end
        end

        # Write attributes and units for dimensions
        Dictu = output_dict()

        for nc in cat(ncs_progn,ncs_diagn,dims=1)
            if nc != 0
                NetCDF.putatt(nc,"global",Dictu)
                NetCDF.putatt(nc,"t",Dict("units"=>"s","long_name"=>"time"))
                NetCDF.putatt(nc,"x",Dict("units"=>"m","long_name"=>"zonal coordinate"))
                NetCDF.putatt(nc,"y",Dict("units"=>"m","long_name"=>"meridional coordinate"))
            end
        end

        # write initial conditions
        iout = 1   # counter for output time steps
        ncs_diagn = output_diagn_nc(ncs_diagn,0,iout,Bu,Bv,LLu1,LLu2,LLv1,LLv2)
        ncs_progn,iout = output_progn_nc(ncs_progn,0,iout,u,v,η,sst)

        # also output scripts
        scripts_output()

        return ncs_progn,ncs_diagn,iout
    else
        return nothing,nothing,nothing
    end
end

function nccreate(x::Array{Float64,1},y::Array{Float64,1},name::String,path::String,unit::String,long_name::String)
    xdim = NcDim("x",length(x),values=x)
    ydim = NcDim("y",length(y),values=y)
    tdim = NcDim("t",0,unlimited=true)

    var = NcVar(name,[xdim,ydim,tdim],t=Float32)
    tvar = NcVar("t",tdim,t=Int32)

    nc = NetCDF.create(path*name*".nc",[var,tvar],mode=NC_NETCDF4)
    NetCDF.putatt(nc,name,Dict("units"=>unit,"long_name"=>long_name))
    return nc
end

"""Writes prognostic variables to pre-initialised netCDF file."""
function output_progn_nc(ncs,i,iout,u,v,η,sst)

    # if nprocs > 1
    #     #TODO MPI Gather data
    #     #TODO rename u,v,η necessary? To distinguish between the loval u,v,η and the gathered u,v,η?
    # end

    # output only every nout time steps
    # only process 0 will do the output
    if i % nout == 0 && output #&& prank == 0

        # cut off the halo
        if ncs[1] != 0
            NetCDF.putvar(ncs[1],"u",Float32.(u[halo+1:end-halo,halo+1:end-halo]),start=[1,1,iout],count=[-1,-1,1])
        end
        if ncs[2] != 0
            NetCDF.putvar(ncs[2],"v",Float32.(v[halo+1:end-halo,halo+1:end-halo]),start=[1,1,iout],count=[-1,-1,1])
        end
        if ncs[3] != 0
            NetCDF.putvar(ncs[3],"eta",Float32.(η[haloη+1:end-haloη,haloη+1:end-haloη]),start=[1,1,iout],count=[-1,-1,1])
        end
        if ncs[4] != 0
            NetCDF.putvar(ncs[4],"sst",Float32.(sst[halosstx+1:end-halosstx,halossty+1:end-halossty]),start=[1,1,iout],count=[-1,-1,1])
        end

        for nc in ncs
            if nc !=0
                #TODO check whether Int64 here clashes with the Int32 of type of time dimension
                NetCDF.putvar(nc,"t",Int64[i*dtint],start=[iout])
                NetCDF.sync(nc)     # sync to view netcdf while model is still running
            end
        end

        iout += 1
    end

    #TODO MPI Barrier, Waitall?

    return ncs,iout
end

""" Writes data to a netCDF file."""
function output_diagn_nc(ncs,i,iout,Bu,Bv,LLu1,LLu2,LLv1,LLv2)
    if i % nout == 0 && output && output_diagn #&& prank == 0

        # cut off the halo
        if ncs[1] != 0
            NetCDF.putvar(ncs[1],"Bu",Float32.(Bu[2-ep:end-1,2:end-1]),start=[1,1,iout],count=[-1,-1,1])
        end
        if ncs[2] != 0
            NetCDF.putvar(ncs[2],"Bv",Float32.(Bv[2:end-1,2:end-1]),start=[1,1,iout],count=[-1,-1,1])
        end
        if ncs[3] != 0
            NetCDF.putvar(ncs[3],"LLu1",Float32.(LLu1[:,2:end-1]),start=[1,1,iout],count=[-1,-1,1])
        end
        if ncs[4] != 0
            NetCDF.putvar(ncs[4],"LLu2",Float32.(LLu2[2-ep:end,:]),start=[1,1,iout],count=[-1,-1,1])
        end
        if ncs[5] != 0
            NetCDF.putvar(ncs[5],"LLv1",Float32.(LLv1[:,2:end-1]),start=[1,1,iout],count=[-1,-1,1])
        end
        if ncs[6] != 0
            NetCDF.putvar(ncs[6],"LLv2",Float32.(LLv2[2:end-1,:]),start=[1,1,iout],count=[-1,-1,1])
        end

        for nc in ncs
            if nc !=0
                #TODO check whether Int64 here clashes with the Int32 of type of time dimension
                NetCDF.putvar(nc,"t",Int64[i*dtint],start=[iout])
                NetCDF.sync(nc)     # sync to view netcdf while model is still running
            end
        end
    end

    #TODO MPI Barrier, Waitall?

    return ncs
end

"""Closes netCDF and progress.txt files."""
function output_close(ncs_progn,ncs_diagn,progrtxt)
    if output #&& prank == 0
        for nc in cat(ncs_progn,ncs_diagn,dims=1)
            if nc !=0
                NetCDF.close(nc)
            end
        end
        println("All data stored.")
        write(progrtxt,"All data stored.")
        close(progrtxt)
    end
end

"""Checks output folders to determine a 4-digit run id number."""
function get_run_id_path(order="continue",run_id=nothing)

    """Finds the first gap in a list of integers."""
    function gap(a::Array{Int,1})
        try
            return minimum([i for i in minimum(a):maximum(a) if ~(i in a)])
        catch
            return maximum(a)+1
        end
    end

    # only process rank 0 checks existing folders
    if output #&& prank == 0
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

#TODO in ensemble mode, the .jl files might have changed since the start and do not correspond to what
#TODO is actually executed!
"""Archives all .jl files of juls in the output folder to make runs reproducible."""
function scripts_output()
    if output #&& prank == 0
        # copy all files in juls main folder
        mkdir(runpath*"scripts")
        for juliafile in filter(x->endswith(x,".jl"),readdir())
            cp(juliafile,runpath*"scripts/"*juliafile)
        end

        # and also in the src folder
        mkdir(runpath*"scripts/src")
        for juliafile in filter(x->endswith(x,".jl"),readdir("src"))
            cp("src/"*juliafile,runpath*"scripts/src/"*juliafile)
        end
    end
end

"""Creates a dictionary with many parameter constants to be included in the nc files."""
function output_dict()
    # Attributes for nc
    Dictu = Dict{String,Any}("description"=>"Data from shallow-water model juls.")
    Dictu["details"] = "Cartesian coordinates, f or beta-plane, Arakawa C-grid"
    Dictu["reference"] = "github.com/milankl/juls"

    Dictu["nx"] = nx
    Dictu["Lx"] = Lx
    Dictu["L_ratio"] = L_ratio
    Dictu["delta"] = Δ

    Dictu["halo"] = halo
    Dictu["haloeta"] = haloη
    Dictu["halosstx"] = halosstx
    Dictu["halossty"] = halossty

    Dictu["g"] = gravity
    Dictu["water_depth"] = water_depth
    Dictu["phi"] = ϕ
    Dictu["density"] = ρ

    Dictu["wind_forcing"] = wind_forcing
    Dictu["Fx0"] = Fx0

    Dictu["topography_feature"] = topography_feature
    Dictu["topofeat_height"] = topofeat_height
    Dictu["topofeat_width"] = topofeat_width

    Dictu["surface_forcing"] = string(surface_forcing)
    Dictu["t_relax"] = t_relax
    Dictu["eta_refh"] = η_refh
    Dictu["η_refw"] = η_refw

    Dictu["Numtype"] = string(Numtype)
    Dictu["output_dt"] = output_dt
    Dictu["nout"] = nout
    Dictu["nadvstep"] = nadvstep
    Dictu["nstep_diff"] = nstep_diff
    Dictu["nstep_advcor"] = nstep_advcor

    Dictu["RKo"] = RKo
    Dictu["cfl"] = cfl
    Dictu["Ndays"] = Ndays

    Dictu["bc_x"] = bc_x
    Dictu["lbc"] = lbc

    Dictu["adv_scheme"] = adv_scheme
    Dictu["dynamics"] = dynamics

    Dictu["bottom_friction"] = bottom_friction
    Dictu["drag"] = drag
    Dictu["taudrag"] = τdrag

    Dictu["diffusion"] = diffusion
    Dictu["nuConst"] = ν_const
    Dictu["c_smag"] = c_smag

    Dictu["tracer_advcetion"] = string(tracer_advection)
    Dictu["tracer_relaxation"] = string(tracer_relaxation)
    Dictu["injection_region"] = injection_region
    Dictu["sstrestart"] = string(sstrestart)
    Dictu["Uadv"] = Uadv
    Dictu["SSTmax"] = SSTmax
    Dictu["SSTmin"] = SSTmin
    Dictu["tauSST"] = τSST
    Dictu["SSTw"] = SSTw
    Dictu["SSTphi"] = SSTϕ

    Dictu["initial_cond"] = initial_cond
    Dictu["init_run_id"] = init_run_id
    Dictu["initpath"] = initpath

    return Dictu
end

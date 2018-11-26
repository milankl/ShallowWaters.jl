"""initialises netCDF files for data output."""
function output_ini(u,v,η,sst)
    # only process with rank 0 defines the netCDF file
    if output == 1 #&& prank == 0

        # Attributes
        Dictu = Dict{String,Any}("description"=>"Data from shallow-water model juls.")
        Dictu["details"] = "Cartesian coordinates, f or beta-plane, Arakawa C-grid"
        Dictu["reference"] = "github.com/milankl/juls"
        Dictu["cfl"] = cfl
        Dictu["g"] = gravity
        Dictu["water_depth"] = water_depth
        Dictu["bc_x"] = bc_x
        Dictu["lbc"] = lbc
        Dictu["drag"] = drag
        Dictu["c_smag"] = c_smag
        Dictu["initial_cond"] = initial_cond
        Dictu["init_run_id"] = init_run_id
        Dictu["phi"] = ϕ
        Dictu["density rho"] = ρ

        Dictu["wind_forcing"] = wind_forcing
        Dictu["Fx0"] = Fx0

        Dictu["surface_forcing"] = string(surface_forcing)
        Dictu["t_relax"] = t_relax
        Dictu["eta_refh"] = η_refh
        Dictu["η_refw"] = η_refw

        Dictu["topography_feature"] = topography_feature
        Dictu["topofeat_height"] = topofeat_height
        Dictu["topofeat_width"] = topofeat_width

        Dictu["Numtype"] = string(Numtype)
        Dictu["output_dt"] = nout*dtint

        if "u" in output_vars
            xudim = NcDim("x",nux,values=x_u)
            yudim = NcDim("y",nuy,values=y_u)
            tdim = NcDim("t",0,unlimited=true)

            uvar = NcVar("u",[xudim,yudim,tdim],t=Float32)
            tvaru = NcVar("t",tdim,t=Int64)

            ncu = NetCDF.create(runpath*"u.nc",[uvar,tvaru],mode=NC_NETCDF4)
            NetCDF.putatt(ncu,"u",Dict("units"=>"m/s","long_name"=>"zonal velocity"))
        else
            ncu = 0
        end

        if "v" in output_vars
            xvdim = NcDim("x",nvx,values=x_v)
            yvdim = NcDim("y",nvy,values=y_v)
            tdim = NcDim("t",0,unlimited=true)

            vvar = NcVar("v",[xvdim,yvdim,tdim],t=Float32)
            tvarv = NcVar("t",tdim,t=Int64)

            ncv = NetCDF.create(runpath*"v.nc",[vvar,tvarv],mode=NC_NETCDF4)
            NetCDF.putatt(ncv,"v",Dict("units"=>"m/s","long_name"=>"meridional velocity"))
        else
            ncv = 0
        end

        if "eta" in output_vars
            xTdim = NcDim("x",nx,values=x_T)
            yTdim = NcDim("y",ny,values=y_T)
            tdim = NcDim("t",0,unlimited=true)

            ηvar = NcVar("eta",[xTdim,yTdim,tdim],t=Float32)
            tvarη = NcVar("t",tdim,t=Int64)

            ncη = NetCDF.create(runpath*"eta.nc",[ηvar,tvarη],mode=NC_NETCDF4)
            NetCDF.putatt(ncη,"eta",Dict("units"=>"m","long_name"=>"sea surface height"))
        else
            ncη = 0
        end

        if "sst" in output_vars
            xTdim = NcDim("x",nx,values=x_T)
            yTdim = NcDim("y",ny,values=y_T)
            tdim = NcDim("t",0,unlimited=true)

            ηvar = NcVar("sst",[xTdim,yTdim,tdim],t=Float32)
            tvarη = NcVar("t",tdim,t=Int64)

            ncsst = NetCDF.create(runpath*"sst.nc",[ηvar,tvarη],mode=NC_NETCDF4)
            NetCDF.putatt(ncsst,"sst",Dict("units"=>"1","long_name"=>"sea surface temperature"))
        else
            ncsst = 0
        end

        ncs = (ncu,ncv,ncη,ncsst)

        # Write attributes and units
        for nc in ncs
            if nc != 0
                NetCDF.putatt(nc,"global",Dictu)
                NetCDF.putatt(nc,"t",Dict("units"=>"s","long_name"=>"time"))
                NetCDF.putatt(nc,"x",Dict("units"=>"m","long_name"=>"zonal coordinate"))
                NetCDF.putatt(nc,"y",Dict("units"=>"m","long_name"=>"meridional coordinate"))
            end
        end

        # write initial conditions
        iout = 1   # counter for output time steps
        ncs,iout = output_nc(ncs,u,v,η,sst,0,iout)

        # also output scripts
        scripts_output()

        return ncs,iout
    else
        return nothing, nothing
    end
end

""" Writes data to a netCDF file."""
function output_nc(ncs,u,v,η,sst,i,iout)

    # if nprocs > 1
    #     #TODO MPI Gather data
    #     #TODO rename u,v,η necessary? To distinguish between the loval u,v,η and the gathered u,v,η?
    # end

    # output only every nout time steps
    # only process 0 will do the output
    if i % nout == 0 && output == 1 #&& prank == 0

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
                NetCDF.putvar(nc,"t",Int64[i*dtint],start=[iout])
                NetCDF.sync(nc)     # sync to view netcdf while model is still running
            end
        end

        iout += 1
    end

    #TODO MPI Barrier, Waitall?

    return ncs,iout
end

"""Closes netCDF and progress.txt files."""
function output_close(ncs,progrtxt)
    if output == 1 #&& prank == 0
        for nc in ncs
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

    function gap(a::Array{Int,1})
        try
            return minimum([i for i in minimum(a):maximum(a) if ~(i in a)])
        catch
            return maximum(a)+1
        end
    end

    # only process rank 0 checks existing folders
    if output == 1 #&& prank == 0
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
    if output == 1 #&& prank == 0
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

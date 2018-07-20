function output_ini(u,v,η)
    if output == 1
        xudim = NcDim("x",nux,values=x_u)
        yudim = NcDim("y",nuy,values=y_u)
        xvdim = NcDim("x",nvx,values=x_v)
        yvdim = NcDim("y",nvy,values=y_v)
        xTdim = NcDim("x",nx,values=x_T)
        yTdim = NcDim("y",ny,values=y_T)
        tdim = NcDim("t",0,unlimited=true)

        uvar = NcVar("u",[xudim,yudim,tdim],t=Float32)
        vvar = NcVar("v",[xvdim,yvdim,tdim],t=Float32)
        ηvar = NcVar("eta",[xTdim,yTdim,tdim],t=Float32)
        tvar = NcVar("t",tdim,t=Int64)

        ncu = NetCDF.create(runpath*"u.nc",[uvar,tvar],mode=NC_NETCDF4)
        ncv = NetCDF.create(runpath*"v.nc",[vvar,tvar],mode=NC_NETCDF4)
        ncη = NetCDF.create(runpath*"eta.nc",[ηvar,tvar],mode=NC_NETCDF4)

        NetCDF.putvar(ncu,"u",Float32.(u),start=[1,1,1],count=[-1,-1,1])
        NetCDF.putvar(ncv,"v",Float32.(v),start=[1,1,1],count=[-1,-1,1])
        NetCDF.putvar(ncη,"eta",Float32.(η),start=[1,1,1],count=[-1,-1,1])

        NetCDF.putvar(ncu,"t",Array(0:nout_total-1)*dtint*nout)
        #NetCDF.putvar(ncv,"t",Array(0:nout_total-1)*dtint*nout)
        #NetCDF.putvar(ncη,"t",Array(0:nout_total-1)*dtint*nout)

        # Attributes
        Dictu = Dict{String,Any}("description"=>"Data from shallow-water model juls.")
        Dictu["details"] = "Cartesian coordinates, f or beta-plane, Arakawa C-grid"
        Dictu["reference"] = "github.com/milankl/juls"
        Dictu["cfl"] = cfl
        Dictu["g"] = gravity
        Dictu["water_depth"] = water_depth
        Dictu["bc_x"] = bc_x
        Dictu["lbc"] = lbc
        Dictu["c_D"] = c_D
        Dictu["c_smag"] = c_smag
        Dictu["initial_cond"] = initial_cond
        Dictu["init_run_id"] = init_run_id
        Dictu["phi"] = ϕ
        Dictu["seamount_height"] = 150.
        Dictu["Numtype"] = string(Numtype)
        Dictu["output_dt"] = nout*dtint

        NetCDF.putatt(ncu,"global",Dictu)
        NetCDF.putatt(ncv,"global",Dictu)
        NetCDF.putatt(ncη,"global",Dictu)

        NetCDF.putatt(ncu,"t",Dict("units"=>"s","long_name"=>"time"))
        NetCDF.putatt(ncv,"t",Dict("units"=>"s","long_name"=>"time"))
        NetCDF.putatt(ncη,"t",Dict("units"=>"s","long_name"=>"time"))

        NetCDF.putatt(ncu,"x",Dict("units"=>"m","long_name"=>"zonal coordinate"))
        NetCDF.putatt(ncv,"x",Dict("units"=>"m","long_name"=>"zonal coordinate"))
        NetCDF.putatt(ncη,"x",Dict("units"=>"m","long_name"=>"zonal coordinate"))

        NetCDF.putatt(ncu,"y",Dict("units"=>"m","long_name"=>"meridional coordinate"))
        NetCDF.putatt(ncv,"y",Dict("units"=>"m","long_name"=>"meridional coordinate"))
        NetCDF.putatt(ncη,"y",Dict("units"=>"m","long_name"=>"meridional coordinate"))

        NetCDF.putatt(ncu,"u",Dict("units"=>"m/s","long_name"=>"zonal velocity"))
        NetCDF.putatt(ncv,"v",Dict("units"=>"m/s","long_name"=>"meridional velocity"))
        NetCDF.putatt(ncη,"eta",Dict("units"=>"m","long_name"=>"sea surface height"))

        iout = 2    # counter for output time steps

        # also output scripts
        scripts_output()

        return (ncu,ncv,ncη),iout
    else
        return nothing, nothing
    end
end

function output_nc(ncs,u,v,η,i,iout)
    if i % nout == 0    # output only every nout time steps
        if output == 1
            NetCDF.putvar(ncs[1],"u",Float32.(u),start=[1,1,iout],count=[-1,-1,1])
            #NetCDF.putvar(ncs[1],"t",[i*dtint],start=[iout],count=[1])
            NetCDF.putvar(ncs[2],"v",Float32.(v),start=[1,1,iout],count=[-1,-1,1])
            NetCDF.putvar(ncs[3],"eta",Float32.(η),start=[1,1,iout],count=[-1,-1,1])
            #println("Time step $iout written to file.")
            iout += 1
        end
    end

    return ncs,iout
end

function output_close(ncs,progrtxt)
    if output == 1
        for nc in ncs
            NetCDF.close(nc)
        end
        println("All data stored.")
        write(progrtxt,"All data stored.")
        close(progrtxt)
    end
end

function get_run_id_path()
        if output == 1
                runlist = filter(x->startswith(x,"run"),readdir(outpath))
                existing_runs = [parse(Int,id[4:end]) for id in runlist]
                if length(existing_runs) == 0           # if no runfolder exists yet
                        runpath = outpath*"run0000/"
                        mkdir(runpath)
                        return 0,runpath
                else
                        run_id = maximum(existing_runs)+1
                        runpath = outpath*"run"*@sprintf("%04d",run_id)*"/"
                        mkdir(runpath)
                        return run_id,runpath
                end
        else
                return 0,"no runpath"
        end
end

const run_id,runpath = get_run_id_path()

function scripts_output()
        if output == 1
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

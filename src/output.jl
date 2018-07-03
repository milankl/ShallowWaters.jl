function output_nc_ini(u,v,η)

    xudim = NcDim("x",nux,values=x_u)
    yudim = NcDim("y",nuy,values=y_u)
    xvdim = NcDim("x",nvx,values=x_v)
    yvdim = NcDim("y",nvy,values=y_v)
    xTdim = NcDim("x",nx,values=x_T)
    yTdim = NcDim("y",ny,values=y_T)
    tdim = NcDim("t",nout_total,values=t_vec,unlimited=true)

    uvar = NcVar("u",[xudim,yudim,tdim],t=Float32)
    vvar = NcVar("v",[xvdim,yvdim,tdim],t=Float32)
    ηvar = NcVar("eta",[xTdim,yTdim,tdim],t=Float32)

    runpath = "data/run0000/"
    ncu = NetCDF.create(runpath*"u.nc",uvar,mode=NC_NETCDF4)
    ncv = NetCDF.create(runpath*"v.nc",vvar,mode=NC_NETCDF4)
    ncη = NetCDF.create(runpath*"eta.nc",ηvar,mode=NC_NETCDF4)

    NetCDF.putvar(ncu,"u",u,start=[1,1,1],count=[-1,-1,1])
    NetCDF.putvar(ncv,"v",v,start=[1,1,1],count=[-1,-1,1])
    NetCDF.putvar(ncη,"eta",η,start=[1,1,1],count=[-1,-1,1])

    iout = 2    # counter for output time steps
    return (ncu,ncv,ncη),iout
end

function output_nc(ncs,u,v,η,i,iout)
    if i % nout == 0    # output only every nout time steps
        if output == 1
            NetCDF.putvar(ncs[1],"u",u,start=[1,1,iout],count=[-1,-1,1])
            NetCDF.putvar(ncs[2],"v",v,start=[1,1,iout],count=[-1,-1,1])
            NetCDF.putvar(ncs[3],"eta",η,start=[1,1,iout],count=[-1,-1,1])
            iout += 1
        end
    end

    return ncs,iout
end

function output_nc_close(ncs)
    if output == 1
        for nc in ncs
            NetCDF.close(nc)
        end
        println("All data stored.")
    end
end

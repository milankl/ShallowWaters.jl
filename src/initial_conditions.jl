function initial_conditions()

    if initial_cond == "rest"
        u = zeros(Numtype,nux,nuy)
        v = zeros(Numtype,nvx,nvy)
        η = zeros(Numtype,nx,ny)

    elseif initial_cond == "ncfile"
        inipath = outpath*"run"*@sprintf("%04d",init_run_id)*"/"

        ncu = NetCDF.open(inipath*"u.nc")
        u = ncu.vars["u"][:,:,end]
        NetCDF.close(ncu)

        ncv = NetCDF.open(inipath*"v.nc")
        v = ncv.vars["v"][:,:,end]
        NetCDF.close(ncv)

        ncη = NetCDF.open(inipath*"eta.nc")
        η = ncη.vars["eta"][:,:,end]
        NetCDF.close(ncη)

        # remove singleton time dimension
        # and convert from Float32 to Numtype
        u = Numtype.(reshape(u,size(u)[1:2]))
        v = Numtype.(reshape(v,size(v)[1:2]))
        η = Numtype.(reshape(η,size(η)[1:2]))
    end

    return u,v,η
end

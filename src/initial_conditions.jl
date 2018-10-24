function initial_conditions()
    # initialise the state matrices u,v,η and set their initial conditions

    if initial_cond == "rest"
        # sub domain sub_dom is the entire domain for a single process (no domain decomposition)
        u = zeros(Numtype,sub_dom["nux"],sub_dom["nuy"])
        v = zeros(Numtype,sub_dom["nvx"],sub_dom["nvy"])
        η = zeros(Numtype,sub_dom["nx"],sub_dom["ny"])

    elseif initial_cond == "ncfile"
        #TODO for domain decomposition: slice ncfile for different processors
        inipath = outpath*"run"*@sprintf("%04d",init_run_id)*"/"

        # take last time step from existing netcdf files
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

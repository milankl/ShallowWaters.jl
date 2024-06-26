struct NcFiles
    u::Union{NcFile,Nothing}        # zonal velocity
    v::Union{NcFile,Nothing}        # meridional velocity
    η::Union{NcFile,Nothing}        # sea surface height anomaly
    sst::Union{NcFile,Nothing}      # tracer / sea surface temperature
    q::Union{NcFile,Nothing}        # potentival vorticity
    ζ::Union{NcFile,Nothing}        # relative vorticity
    du::Union{NcFile,Nothing}       # tendency of u [m^2/s^2]
    dv::Union{NcFile,Nothing}       # tendency of v [m^2/s^2]
    dη::Union{NcFile,Nothing}       # tendency of η [m^2/s]
    H::Union{NcFile,Nothing}        # layer thickness at rest [m]
    iout::Array{Integer,1}          # output index, stored in array for mutability
end

"""Generator function for "empty" NcFiles struct."""
NcFiles(x::Nothing) = NcFiles(x,x,x,x,x,x,x,x,x,x,[0])

"""Generator function for NcFiles struct, creating the underlying netCDF files."""
function NcFiles(feedback::Feedback,S::ModelSetup)

    if S.parameters.output

        @unpack output_vars,compression_level = S.parameters
        P = S.parameters
        @unpack x_u,x_v,x_T,x_q = S.grid
        @unpack y_u,y_v,y_T,y_q = S.grid

        @unpack run_id,runpath = feedback

        ncu = if "u" in output_vars nc_create(x_u,y_u,"u",runpath,"m/s","zonal velocity",P) else nothing end
        ncv = if "v" in output_vars nc_create(x_v,y_v,"v",runpath,"m/s","meridional velocity",P) else nothing end
        ncη = if "η" in output_vars nc_create(x_T,y_T,"eta",runpath,"m","sea surface height",P) else nothing end
        ncsst = if "sst" in output_vars nc_create(x_T,y_T,"sst",runpath,"1","sea surface temperature",P) else nothing end
        ncq = if "q" in output_vars nc_create(x_q,y_q,"q",runpath,"1/(ms)","potential vorticity",P) else nothing end
        ncζ = if "ζ" in output_vars nc_create(x_q,y_q,"relvort",runpath,"1","relative vorticity",P) else nothing end
        ncdu = if "du" in output_vars nc_create(x_u,y_u,"du",runpath,"m^2/s^2","zonal velocity tendency",P) else nothing end
        ncdv = if "dv" in output_vars nc_create(x_v,y_v,"dv",runpath,"m^2/s^2","meridional velocity tendency",P) else nothing end
        ncdη = if "dη" in output_vars nc_create(x_T,y_T,"deta",runpath,"m^2/s","sea surface height tendency",P) else nothing end
        nH = if "H" in output_vars nc_create(x_T,y_T,"H",runpath,"m","undisturbed layer thickness",P,false) else nothing end

        for nc in (ncu,ncv,ncη,ncsst,ncq,ncζ,ncdu,ncdv,ncdη)
            if !isnothing(nc)
                NetCDF.putatt(nc,"t",Dict("units"=>"s","long_name"=>"time"))
                NetCDF.putatt(nc,"x",Dict("units"=>"m","long_name"=>"zonal coordinate"))
                NetCDF.putatt(nc,"y",Dict("units"=>"m","long_name"=>"meridional coordinate"))
            end
        end

        if !isnothing(nH) # output constant fields 
            NetCDF.putatt(nH,"x",Dict("units"=>"m","long_name"=>"zonal coordinate"))
            NetCDF.putatt(nH,"y",Dict("units"=>"m","long_name"=>"meridional coordinate"))
        end

        return NcFiles(ncu,ncv,ncη,ncsst,ncq,ncζ,ncdu,ncdv,ncdη,nH,[0])
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
                    long_name::String,
                    P::Parameter,
                    output_tdim::Bool=true) where {T<:Real}

    @unpack compression_level = P

    xdim = NcDim("x",length(x),values=x)
    ydim = NcDim("y",length(y),values=y)
    
    if output_tdim 
        tdim = NcDim("t",0,unlimited=true)
        var = NcVar(name,[xdim,ydim,tdim],t=Float32,compress=compression_level)
        tvar = NcVar("t",tdim,t=Int32)
        nc = NetCDF.create(joinpath(path,name*".nc"),[var,tvar],mode=NC_NETCDF4)
    else 
        var = NcVar(name,[xdim,ydim],t=Float32,compress=compression_level)
        nc = NetCDF.create(joinpath(path,name*".nc"),[var],mode=NC_NETCDF4)
    end 

    # add missing_value although irrelevant for ncview compatibility
    NetCDF.putatt(nc,name,Dict("units"=>unit,"long_name"=>long_name,"missing_value"=>-999999f0))
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
        @unpack f_q,ep,dtint = S.grid
        @unpack scale,scale_inv,scale_sst = S.constants

        # CUT OFF HALOS
        # As output is before copyto!(u,u0), take u0,v0,η0
        # Tendencies calculate from the last time step, du = u_n+1-u_n etc
        # WRITING THE VARIABLES
        if !isnothing(ncs.u)
            @views u = Float32.(scale_inv*Diag.RungeKutta.u0[halo+1:end-halo,halo+1:end-halo])
            NetCDF.putvar(ncs.u,"u",u,start=[1,1,iout],count=[-1,-1,1])
        end
        if !isnothing(ncs.v)
            @views v = Float32.(scale_inv*Diag.RungeKutta.v0[halo+1:end-halo,halo+1:end-halo])
            NetCDF.putvar(ncs.v,"v",v,start=[1,1,iout],count=[-1,-1,1])
        end
        if !isnothing(ncs.η)
            @views η = Float32.(Diag.RungeKutta.η0[haloη+1:end-haloη,haloη+1:end-haloη])
            NetCDF.putvar(ncs.η,"eta",η,start=[1,1,iout],count=[-1,-1,1])
        end
        if !isnothing(ncs.sst)
            @views sst = Float32.(Prog.sst[halosstx+1:end-halosstx,halossty+1:end-halossty]/scale_sst)
            NetCDF.putvar(ncs.sst,"sst",sst,start=[1,1,iout],count=[-1,-1,1])
        end
        if !isnothing(ncs.q)
            @views q = Float32.(scale_inv*Diag.Vorticity.q[haloη+1:end-haloη,haloη+1:end-haloη])
            NetCDF.putvar(ncs.q,"q",q,start=[1,1,iout],count=[-1,-1,1])
        end
        if !isnothing(ncs.ζ)
            @unpack dvdx,dudy = Diag.Vorticity
            @unpack f_q = S.grid
            @views ζ = Float32.((dvdx[2:end-1,2:end-1]-dudy[2+ep:end-1,2:end-1])./abs.(f_q))
            NetCDF.putvar(ncs.ζ,"relvort",ζ,start=[1,1,iout],count=[-1,-1,1])
        end
        if !isnothing(ncs.du)
            @views u = Float32.(scale_inv*Diag.RungeKutta.u0[halo+1:end-halo,halo+1:end-halo])
            @views du = u-Float32.(scale_inv*Prog.u[halo+1:end-halo,halo+1:end-halo])
            NetCDF.putvar(ncs.du,"du",du,start=[1,1,iout],count=[-1,-1,1])
        end
        if !isnothing(ncs.dv)
            @views v = Float32.(scale_inv*Diag.RungeKutta.v0[halo+1:end-halo,halo+1:end-halo])
            @views dv = v-Float32.(scale_inv*Prog.v[halo+1:end-halo,halo+1:end-halo])
            NetCDF.putvar(ncs.dv,"dv",dv,start=[1,1,iout],count=[-1,-1,1])
        end
        if !isnothing(ncs.dη)
            @views η = Float32.(Diag.RungeKutta.η0[haloη+1:end-haloη,haloη+1:end-haloη])
            @views dη = η-Float32.(Prog.η[haloη+1:end-haloη,haloη+1:end-haloη])
            NetCDF.putvar(ncs.dη,"deta",dη,start=[1,1,iout],count=[-1,-1,1])
        end
        if i==0 && !isnothing(ncs.H) # constant field, output only once, before the timestepping
            @views H = Float32.(S.forcing.H[haloη+1:end-haloη,haloη+1:end-haloη])
            NetCDF.putvar(ncs.H,"H",H,start=[1,1],count=[-1,-1])
        end

        # WRITING THE TIME
        for nc in (ncs.u,ncs.v,ncs.η,ncs.sst,ncs.q,ncs.ζ,ncs.du,ncs.dv,ncs.dη)
            if !isnothing(nc)
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
function get_run_id_path(S::ModelSetup)

    @unpack output,outpath,get_id_mode = S.parameters

    if output
        
        pattern = r"run_\d\d\d\d"               # run_???? in regex
        runlist = filter(x->startswith(x,pattern),readdir(outpath))
        runlist = filter(x->endswith(  x,pattern),runlist)
        existing_runs = [parse(Int,id[5:end]) for id in runlist]

        # get the run id from existing folders
        if length(existing_runs) == 0           # if no runfolder exists yet
            run_id = 1                          # start with run_0001
        else
            run_id = maximum(existing_runs)+1   # next run gets id +1
        end
        
        id = @sprintf("%04d",run_id)

        run_id2 = string("run_",id)
        run_path = joinpath(outpath,run_id2)
        @assert !(run_id2 in readdir(outpath)) "Run folder $run_path already exists."
        mkdir(run_path)             # actually create the folder

        return run_id, run_path

    else
        return 0,"no runpath"
    end
end
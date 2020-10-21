using NetCDF

x = Array{Float32}(1:10)
y = Array{Float32}(1:10)

xdim = NcDim("x",length(x),values=x)
ydim = NcDim("y",length(y),values=y)
tdim = NcDim("t",0,unlimited=true)

var = NcVar("u",[xdim,ydim,tdim],t=Float32)
tvar = NcVar("t",tdim,t=Int32)

nc = NetCDF.create("uu.nc",[var,tvar],mode=NC_NETCDF4)
NetCDF.putatt(nc,"u",Dict("units"=>"m/s","long_name"=>"velocity"))

for i in 1:10
    NetCDF.putvar(nc,"u",rand(Float32,length(x),length(y)),start=[1,1,i],count=[-1,-1,1])
end

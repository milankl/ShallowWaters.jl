include("Juls.jl")
using .Juls
using BFloat16s

Base.round(x::BFloat16, r::RoundingMode{:Up}) = BFloat16(ceil(Float32(x)))
Base.round(x::BFloat16, r::RoundingMode{:Down}) = BFloat16(floor(Float32(x)))
Base.round(x::BFloat16, r::RoundingMode{:Nearest}) = BFloat16(round(Float32(x)))

a = BFloat16(1.2)
floor(a)

#RunJuls(BFloat16,Tprog=Float32)
#RunJuls(BFloat16,Tprog=Float32,output=true,nx=200,Ndays=100,initial_cond="ncfile",sst_initial="rect",sst_rect_coords=[0.,0.5,0.,1.0],outpath="/local/home/kloewer/julsdata/mixedprec")

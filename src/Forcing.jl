struct Forcing{T<:AbstractFloat}
    Fx::Array{T,2}
    Fy::Array{T,2}
    H::Array{T,2}
    η_ref::Array{T,2}
    Fη::Array{T,2}
    sst_ref::Array{T,2}
    sst_γ::Array{T,2}
end

function Forcing{T}(P::Parameter,G::Grid)

    if wind_forcing_x == "channel"
        Fx,_ = ChannelWind(T,P,G)
    elseif wind_forcing_x == "shear"
        Fx,_ = ShearWind(T,P,G)
    elseif wind_forcing_x == "double_gyre"
        Fx,_ = DoubleGyreWind(T,P,G)
    elseif wind_forcing_x == "constant"
        Fx,_ = ConstantWind(T,P,G)
    end

    if wind_forcing_y == "channel"
        _,Fy = ChannelWind(T,P,G)
    elseif wind_forcing_y == "shear"
        _,Fy = ShearWind(T,P,G)
    elseif wind_forcing_y == "double_gyre"
        _,Fy = DoubleGyreWind(T,P,G)
    elseif wind_forcing_y == "constant"
        _,Fy = ConstantWind(T,P,G)
    end



"""Returns the constant forcing matrices Fx,Fy that vary only meriodionally/zonally
as a cosine with strongest forcing in the middle and vanishing forcing at boundaries."""
function ChannelWind(::Type{T},P::Parameter,G::Grid) where {T<:AbstractFloat}

    @unpack Δ,x_u,y_u,Lx,Ly = G
    @unpack Fx0,Fy0,H,ρ = P

    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    xx_u,yy_u = meshgrid(x_u,y_u)
    Fx = (Δ*Fx0/ρ/H)*cos.(π*(yy_u/Ly .- 1/2)).^2
    Fy = (Δ*Fy0/ρ/H)*cos.(π*(xx_u/Lx .- 1/2)).^2

    return T.(Fx),T.(Fy)
end

"""Returns the constant forcing matrices Fx,Fy that vary only meriodionally/zonally
as a hyperbolic tangent with strongest shear in the middle."""
function ShearWind(::Type{T},P::Parameter,G::Grid) where {T<:AbstractFloat}

    @unpack Δ,x_u,y_u,Lx,Ly = G
    @unpack Fx0,Fy0,H,ρ = P

    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    xx_u,yy_u = meshgrid(x_u,y_u)
    Fx = (Δ*Fx0/ρ/H)*tanh.(2π*(yy_u/Ly .- 1/2))
    Fy = (Δ*Fy0/ρ/H)*tanh.(2π*(xx_u/Lx .- 1/2))

    return T.(Fx),T.(Fy)
end

"""Returns the constant forcing matrices Fx,Fy that vary only meriodionally/zonally
with a superposition of sin & cos for a double gyre circulation.
See Cooper&Zanna 2015 or Kloewer et al 2018."""
function DoubleGyreWind(::Type{T},P::Parameter,G::Grid) where {T<:AbstractFloat}

    @unpack Δ,x_u,y_u,Lx,Ly = G
    @unpack Fx0,Fy0,H,ρ = P

    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    xx_u,yy_u = meshgrid(x_u,y_u)
    Fx = (Δ*Fx0/ρ/H)*(cos.(2π*(yy_u/Ly .- 1/2)) + 2*sin.(π*(yy_u/Ly .- 1/2)))
    Fy = (Δ*Fy0/ρ/H)*(cos.(2π*(xx_u/Lx .- 1/2)) + 2*sin.(π*(xx_u/Lx .- 1/2)))
    return T.(Fx),T.(Fy)
end

"""Returns constant in in space forcing matrices Fx,Fy."""
function ConstantWind(::Type{T},P::Parameter,G::Grid) where {T<:AbstractFloat}

    @unpack Δ,x_u,y_u,Lx,Ly = G
    @unpack Fx0,Fy0,H,ρ = P

    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    Fx = (Δ*Fx0/ρ/H)*ones(T,nux,nuy)
    Fy = (Δ*Fy0/ρ/H)*ones(T,nvx,nvy)

    return Fx,Fy
end


"""Returns a reference state for Newtonian cooling/surface relaxation shaped as a
hyperbolic tangent to force the continuity equation."""
function interface_relaxation()
    # width of a tangent is defined as the distance between 65% of its minimum value and 65% of the max
    xx_T,yy_T = meshgrid(x_T,y_T)
    η_ref = -(η_refh/2)*tanh.(2π*(Ly/(4*η_refw))*(yy_T/Ly .- 1/2))
    return Numtype.(η_ref)
end

"""Returns Kelvin wave pumping forcing of the continuity equation at the equator."""
function kelvin_pump(x::AbstractVector,y::AbstractVector)
    β = β_at_lat(ϕ)
    c = wave_speed(gravity,water_depth)
    xx,yy = meshgrid(x,y)


    # y-coordinate of the Equator
    mϕ = m_per_lat()
    y_eq = Ly/2 - ϕ*mϕ
    y_15S = Ly/2 - (ϕ+15)*mϕ
    y_15N = Ly/2 - (ϕ-15)*mϕ

    w_win = 300e3
    Lx_win = 0.8

    Fη = A₀*Δ*exp.(-β*(yy.-y_eq).^2/(2c))#.*
        #(1/2 .+ 1/2*tanh.(2π*(Lx/(4*w_win))*(xx/Lx .- (1-Lx_win)/2))).*
        #(1/2 .- 1/2*tanh.(2π*(Lx/(4*w_win))*(xx/Lx .- (1-(1-Lx_win)/2))))

    Fη[yy .< y_15S] .= 0.0
    Fη[yy .> y_15N] .= 0.0

    return Numtype.(Fη)
end

function Fηt(t::Real)
    return -1*Numtype(sin(2*π*ωyr*t/3600/24/365))
end


"""Returns a matrix of water depth for the whole domain that contains a
Gaussian seamount in the middle. Water depth, heigth and width of the
seamount are adjusted with the constants water_depth, topofeat_height and topofeat_width."""
function seamount()
    xx_T,yy_T = meshgrid(x_T_halo,y_T_halo)
    bumpx = exp.(-((xx_T .- Lx/2).^2)/(2*topofeat_width^2))
    bumpy = exp.(-((yy_T .- Ly/2).^2)/(2*topofeat_width^2))

    H = water_depth .- topofeat_height*bumpx.*bumpy
    return Numtype.(H)
end

"""Returns a matrix of water depth for the whole domain that contains a
meridional Gaussian ridge in the middle. Water depth, heigth and width of the
ridge are adjusted with the constants water_depth, topofeat_height and topofeat_width."""
function ridge()
    xx_T,yy_T = meshgrid(x_T_halo,y_T_halo)
    bumpx = exp.(-((xx_T .- Lx/2).^2)/(2*topofeat_width^2))

    H = water_depth .- topofeat_height*bumpx
    return Numtype.(H)
end

"""Same as ridge() but for 3 ridges at 1/4,1/2,3/4 of the domain."""
function ridges()
    xx_T,yy_T = meshgrid(x_T_halo,y_T_halo)

    # bumps in x direction
    # shift slightly left/right to avoid a symmetric solution
    b0x = exp.(-(xx_T.^2)/(2*topofeat_width^2))
    b1x = exp.(-((xx_T .- 0.99*Lx/4).^2)/(2*topofeat_width^2))
    b2x = exp.(-((xx_T .- 1.01*Lx/2).^2)/(2*topofeat_width^2))
    b3x = exp.(-((xx_T .- 0.99*3*Lx/4).^2)/(2*topofeat_width^2))
    b4x = exp.(-((xx_T .- Lx).^2)/(2*topofeat_width^2))

    th = topofeat_height    # for convenience
    H = water_depth .- th*b0x .- th*b1x .- th*b2x .- th*b3x .- th*b4x
    return Numtype.(H)
end

"""Bathtub"""
function bathtub()

    x = range(-2.0,stop=1.1,length=length(x_T_halo))
    y = range(-1.1,stop=1.1,length=length(y_T_halo))
    xx,yy = meshgrid(x,y)

    B = xx.^10 + yy.^10
    B[B .> 1.0] .= 1.0

    H = water_depth .- topofeat_height*B
    return Numtype.(H)
end

"""Returns a matrix of constant water depth specified by the constant water_depth."""
function flat_bottom()
    H = fill(water_depth,(nx+2*haloη,ny+2*haloη))
    return Numtype.(H)
end

# rename for convenience
if topography_feature == "ridge"
    topography = ridge
elseif topography_feature == "ridges"
    topography = ridges
elseif topography_feature == "seamount"
    topography = seamount
elseif topography_feature == "flat"
    topography = flat_bottom
elseif topography_feature == "bathtub"
    topography = bathtub
else
    throw(error("Topography feature not correctly declared."))
end

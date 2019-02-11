"""Returns the constant forcing matrix Fx that varies only meriodionally
as a cosine with strongest forcing in the middle and vanishing forcing at the northern
and southern boundary."""
function channel_wind()
    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    xx_u,yy_u = meshgrid(x_u,y_u)
    Fx = (Δ*Fx0/ρ/water_depth)*cos.(π*(yy_u/Ly .- 1/2)).^2
    return Numtype.(Fx)
end

"""Returns the constant forcing matrix Fx that varies only meriodionally
as a hyperbolic tangent with strongest shear in the middle."""
function shear_wind()
    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    xx_u,yy_u = meshgrid(x_u,y_u)
    Fx = (Δ*Fx0/ρ/water_depth)*tanh.(2π*(yy_u/Ly .- 1/2))
    return Numtype.(Fx)
end

"""Returns the constant in space forcing matrix Fx."""
function constant_wind()
    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    xx_u,yy_u = meshgrid(x_u,y_u)
    Fx = (Δ*Fx0/ρ/water_depth)*ones(size(xx_u))
    return Numtype.(Fx)
end

"""Returns the constant in space forcing matrix Fy."""
function constant_wind_y()
    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    xx_v,yy_v = meshgrid(x_v,y_v)
    Fy = (Δ*Fy0/ρ/water_depth)*ones(size(xx_v))
    return Numtype.(Fy)
end

"""Returns the constant forcing matrix Fx that varies only meriodionally
with a superposition of sin & cos for a double gyre circulation.
See Cooper&Zanna 2015 or Kloewer et al 2018."""
function double_gyre_wind()
    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    xx_u,yy_u = meshgrid(x_u,y_u)
    Fx = (Δ*Fx0/ρ/water_depth)*(cos.(2π*(yy_u/Ly .- 1/2)) + 2*sin.(π*(yy_u/Ly .- 1/2)))
    return Numtype.(Fx)
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

if wind_forcing_x == "channel"
    windx = channel_wind
elseif wind_forcing_x == "shear"
    windx = shear_wind
elseif wind_forcing_x == "double_gyre"
    windx = double_gyre_wind
elseif wind_forcing_x == "constant"
    windx = constant_wind
else
    throw(error("Wind forcing not correctly specified."))
end

if wind_forcing_y == "constant"
    windy = constant_wind_y
else
    throw(error("Wind forcing not correctly specified."))
end

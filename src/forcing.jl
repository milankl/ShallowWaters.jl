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

"""Returns the constant forcing matrix Fx that varies only meriodionally
with a superposition of sin & cos for a double gyre circulation.
See Cooper&Zanna 2015 or Kloewer et al 2018."""
function double_gyre_wind()
    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    xx_u,yy_u = meshgrid(x_u,y_u)
    Fx = (Δ*Fx0/ρ/water_depth)*(cos.(2π*(yy_u/Ly .- 1/2)) + 2*sin.(π*(yy_u/Ly .- 1/2)))
    return Numtype.(Fx)
end

"""Returns a zero matrix for no wind forcing."""
function no_wind()
    return zeros(Numtype,nux,nuy)
end


"""Returns a reference state for Newtonian cooling/surface relaxation shaped as a
hyperbolic tangent to force the continuity equation."""
function interface_relaxation()
    # width of a tangent is defined as the distance between 65% of its minimum value and 65% of the max
    xx_T,yy_T = meshgrid(x_T,y_T)
    η_ref = -(η_refh/2)*tanh.(2π*(Ly/(4*η_refw))*(yy_T/Ly .- 1/2))
    return Numtype.(η_ref)
end

if wind_forcing == "channel"
    wind = channel_wind
elseif wind_forcing == "shear"
    wind = shear_wind
elseif wind_forcing == "double_gyre"
    wind = double_gyre_wind
elseif wind_forcing == "constant"
    wind = constant_wind
elseif wind_forcing == "none"
    wind = no_wind
else
    throw(error("Wind forcing not correctly specified."))
end

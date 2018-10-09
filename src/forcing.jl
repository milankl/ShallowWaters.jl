function channel_wind()

    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    xx_u,yy_u = meshgrid(x_u,y_u)
    Fx = (Δ*Fx0/ρ/water_depth)*cos.(π*(yy_u/Ly .- 1/2)).^2
    return Numtype.(Fx)
end

function shear_wind()

    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    xx_u,yy_u = meshgrid(x_u,y_u)
    Fx = (Δ*Fx0/ρ/water_depth)*tanh.(2π*(yy_u/Ly .- 1/2))
    return Numtype.(Fx)
end


function double_gyre_wind()

    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    xx_u,yy_u = meshgrid(x_u,y_u)
    Fx = (Δ*Fx0/ρ/water_depth)*(cos.(2π*(yy_u/Ly .- 1/2)) + 2*sin.(π*(yy_u/Ly .- 1/2)))
    return Numtype.(Fx)
end

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
else
    throw(error("Wind forcing not correctly specified. Allowed: 'channel', 'shear', 'double_gyre'"))
end

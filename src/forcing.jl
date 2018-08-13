function channel_wind()

    # amplitude
    Fx0 = 0.12
    xx_u,yy_u = meshgrid(x_u,y_u)

    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    Fx = (Δ*Fx0/ρ/water_depth)*cos.(π*(yy_u/Ly-1/2)).^2
    return Numtype.(Fx)
end

function shear_wind()
    # amplitude
    Fx0 = 0.12
    xx_u,yy_u = meshgrid(x_u,y_u)

    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    Fx = (Δ*Fx0/ρ/water_depth)*tanh.(2π*(yy_u/Ly-1/2))
    return Numtype.(Fx)
end


function double_gyre_wind()

    # amplitude
    Fx0 = 0.12
    xx_u,yy_u = meshgrid(x_u,y_u)

    # for non-dimensional gradients the wind forcing needs to contain the grid spacing Δ
    Fx = (Δ*Fx0/ρ/water_depth)*(cos.(2π*(yy_u/Ly-1/2)) + 2*sin.(π*(yy_u/Ly - 1/2)))
    return Numtype.(Fx)
end

# change wind forcing here
wind = channel_wind
# wind = double_gyre_wind

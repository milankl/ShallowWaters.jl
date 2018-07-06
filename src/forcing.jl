function channel_wind()

    # amplitude
    Fx0 = 0.12
    xx_u,yy_u = meshgrid(x_u,y_u)

    Fx = (Δ*Fx0/ρ/water_depth)*cos.(π*(yy_u/Ly-1/2)).^2
    return Numtype.(Fx)
end


function double_gyre_wind()

    # amplitude
    Fx0 = 0.12
    xx_u,yy_u = meshgrid(x_u,y_u)

    Fx = (Δ*Fx0/ρ/water_depth)*(cos.(2π*(yy_u/Ly-1/2)) + 2*sin.(π*(yy_u/Ly - 1/2)))
    return Numtype.(Fx)
end

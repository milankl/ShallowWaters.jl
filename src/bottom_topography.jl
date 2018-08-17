function seamount()
    # A gaussian seamount
    xx_T,yy_T = meshgrid(x_T_halo,y_T_halo)

    σx = 200e3  # extent in x [m]
    σy = 200e3  # extent in y [m]

    bumpx = exp.(-((xx_T - Lx/2).^2)/(2*σx^2))
    bumpy = exp.(-((yy_T - Ly/2).^2)/(2*σy^2))

    H = water_depth - seamount_height*bumpx.*bumpy
    return Numtype.(H)
end

function ridge()
    # A gaussian ridge
    xx_T,yy_T = meshgrid(x_T_halo,y_T_halo)

    σx = 200e3  # extent in x [m]

    bumpx = exp.(-((xx_T - Lx/2).^2)/(2*σx^2))

    H = water_depth - seamount_height*bumpx
    return Numtype.(H)
end

function flat_bottom()
    # A flat bottom
    xx_T,yy_T = meshgrid(x_T_halo,y_T_halo)
    H = water_depth*ones(xx_T)
    return Numtype.(H)
end

# set the desired topography here
topography = ridge
# topography = seamount
# topography = flat_bottom

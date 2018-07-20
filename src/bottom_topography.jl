function seamount()
    # A gaussian seamount
    xx_T,yy_T = meshgrid(x_T,y_T)

    σx = 200e3
    σy = 200e3

    bumpx = exp.(-((xx_T - Lx/2).^2)/(2*σx^2))
    bumpy = exp.(-((yy_T - Ly/2).^2)/(2*σy^2))

    H = water_depth - seamount_height*bumpx.*bumpy
    return Numtype.(H)
end

function flat_bottom()
    return Numtype(water_depth)
end

function seamount()
    # A gaussian seamount
    xx_T,yy_T = meshgrid(x_T_halo,y_T_halo)

    σx = 200e3  # extent in x [m]
    σy = 200e3  # extent in y [m]

    bumpx = exp.(-((xx_T .- Lx/2).^2)/(2*σx^2))
    bumpy = exp.(-((yy_T .- Ly/2).^2)/(2*σy^2))

    H = water_depth .- seamount_height*bumpx.*bumpy
    return Numtype.(H)
end

function ridge()
    # A gaussian ridge
    xx_T,yy_T = meshgrid(x_T_halo,y_T_halo)

    σx = 200e3  # extent in x [m]

    bumpx = exp.(-((xx_T .- Lx/2).^2)/(2*σx^2))

    H = water_depth .- seamount_height*bumpx
    return Numtype.(H)
end

function flat_bottom()
    # A flat bottom
    H = fill(water_depth,(nx+2*haloη,ny+2*haloη))
    return Numtype.(H)
end

# rename for convenience
if topography_feature == "ridge"
    topography = ridge
elseif topography_feature == "seamount"
    topography = seamount
elseif topography_feature == "flat_bottom"
    topography = flat_bottom
end

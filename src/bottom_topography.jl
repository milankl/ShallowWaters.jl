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

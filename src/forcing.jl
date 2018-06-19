function forcing():
    return 1
end


#=


## FORCING FIELDS
def set_forcing():
    """ Sets up the forcing field Fx of shape Nu and makes them globally available.
    This forcing is constant in time and x-direction, zero in y-direction, and varies only with y.
    Resembles trade winds and westerlies for a double gyre set up taken from Cooper and Zanna, 2015, Ocean Modelling.
    Note: the vector Fx excludes the 1/h-factor which is included in each model time step.
    """
    global Fx

    Lx,Ly = param['Lx'],param['Ly']     # for convenience
    param['rho'] = 1e3                  # density

    xx_u,yy_u = np.meshgrid(param['x_u'],param['y_u'])

    param['Fx0'] = 0.12      # was 0.12
    Fx = param['Fx0']*(np.cos(2*np.pi*(yy_u-Ly/2)/Ly) + 2*np.sin(np.pi*(yy_u - Ly/2)/Ly)) / param['rho']

    # from matrix to vector and set data type to have single/double precision
    Fx = Fx.flatten().astype(param['dat_type'])
=#

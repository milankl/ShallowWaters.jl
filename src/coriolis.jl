function coriolis(lat_0::Float,Ly::Float)

    omega = 2π/(24.*3600.)      # Earth's angular frequency [s^-1]
    R = 6.371e6                 # Earth's radius [m]

    f_0 = 2*omega*sind(lat_0)
    beta = 2*omega/R*cosd(lat_0)

    # subtract the regions mid-y so that phi_0 corresponds to a central latitude
    yy_u = np.array([y_u - Ly/2.]*(nx-1)).T
    yy_v = np.array([y_v - Ly/2.]*nx).T
    yy_q = np.array([y_q - Ly/2.]*(nx+1)).T
    yy_T = np.array([y_T - Ly/2.]*nx).T

    # globally available coriolis parameters (only f_q is actually needed though)
    f_u = (f_0 + beta*yy_u.flatten()).astype(param['dat_type'])
    f_v = (f_0 + beta*yy_v.flatten()).astype(param['dat_type'])
    f_q = (f_0 + beta*yy_q.flatten()).astype(param['dat_type'])
    f_T = (f_0 + beta*yy_T.flatten()).astype(param['dat_type'])

end

function timestep()
    # shallow water gravity wave phase speed
    c_phase = sqrt(g*H)

    # the model timestep dt based on cfl stability criterion to resolve gravity waves
    # converting to integer, i.e. rounding up (ceil)
    dt = cfl*dx/c_phase
    Nt = Int(ceil(Ndays*3600*24/dt)) # number of time steps to integrate
    return dt,Nt
end

function timestep()
    # shallow water gravity wave phase speed
    c_phase = sqrt(g*H)

    # the model timestep dt based on cfl stability criterion to resolve gravity waves
    # converting to integer, i.e. rounding up (ceil)
    dt = cfl*dx/c_phase
    nt = Int(ceil(Ndays*3600*24/dt)) # number of time steps to integrate
    dt = Numtype(dt)
    return dt,nt
end

function time_integration(u,v,η)

    # FORCING
    Fx = double_gyre_wind()
    f_u,f_v,f_q = beta_plane()

    # PREALLOCATE
    du,u0,u1,v_u,dηdx,dLu = preallocate_u_vars(u)
    dv,v0,v1,u_v,dηdy,dLv = preallocate_v_vars(v)
    dη,η0,η1,dudx,dvdy = preallocate_T_variables(η)

    RKa = Numtype.([1/6,1/3,1/3,1/6])
    RKb = Numtype.([.5,.5,1.])

    dt,nt = timestep()

    for i = 1:nt

        u1[:],v1[:],η1[:] = u,v,η

        for rki = 1:4
            du[:],dv[:],dη[:] = rhs(du,dv,dη,u1,v1,η1,Fx,f_u,f_v,v_u,u_v,dηdx,dηdy,dLu,dLv,dudx,dvdy)

            if rki < 4
                u1[:] = u + RKb[rki]*dt*du
                v1[:] = v + RKb[rki]*dt*dv
                η1[:] = η + RKb[rki]*dt*dη
            end

            u0 += RKa[rki]*dt*du
            v0 += RKa[rki]*dt*dv
            η0 += RKa[rki]*dt*dη
        end

        u[:],v[:],η[:] = u0,v0,η0
    end

    return u,v,η
end

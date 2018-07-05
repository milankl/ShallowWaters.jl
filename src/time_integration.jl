function time_integration(u,v,η)

    # FORCING
    Fx = double_gyre_wind()
    f_u,f_v,f_q = beta_plane()

    # PREALLOCATE
    du,u0,u1,dpdx,U,V_u,h_u,q_u,adv_u,dLu,dLu2,νSmag_u = preallocate_u_vars(u)
    dv,v0,v1,dpdy,V,U_v,h_v,q_v,adv_v,dLv,dLv2,νSmag_v = preallocate_v_vars(v)
    dη,η0,η1,dudx,dvdy,h,KEu,KEv,p,νSmag = preallocate_T_variables(η)
    q,h_q,dvdx,dudy,shear = preallocate_q_variables()

    RKa = Numtype.([1/6,1/3,1/3,1/6])
    RKb = Numtype.([.5,.5,1.])

    # feedback and output
    t0 = feedback_ini()
    ncs, iout = output_nc_ini(u,v,η)
    nans_detected = false

    t = 0           # model time
    for i = 1:nt

        u1[:],v1[:],η1[:] = u,v,η

        for rki = 1:4
            rhs(du,dv,dη,u1,v1,η1,Fx,f_q,
                dpdx,dpdy,dLu,dLu2,dLv,dLv2,dudx,dvdy,
                p,KEu,KEv,
                h,h_u,h_v,U,V,U_v,V_u,
                q,dvdx,dudy,h_q,q_u,q_v,
                adv_u,adv_v,
                shear,νSmag,νSmag_u,νSmag_v)

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
        t += dtint

        # feedback and output
        t0,nans_detected = feedback(u,v,η,i,t0,nt,nans_detected)
        ncs,iout = output_nc(ncs,u,v,η,i,iout)
    end

    # feeback and output
    feedback_end(t0)
    output_nc_close(ncs)

    return u,v,η
end

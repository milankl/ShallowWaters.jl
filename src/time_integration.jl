function time_integration(u,v,η)

    # FORCING
    Fx = double_gyre_wind()
    f_u,f_v,f_q = beta_plane()

    # PREALLOCATE
    du,u0,u1,dpdx,U,V_u,h_u,q_u,adv_u,Lu1,Lu2 = preallocate_u_vars(u)
    dv,v0,v1,dpdy,V,U_v,h_v,q_v,adv_v,Lv1,Lv2 = preallocate_v_vars(v)
    dη,η0,η1,dudx,dvdy,h,KEu,KEv,p,νSmag,dLudx,dLvdy,shear = preallocate_T_variables(η)
    q,h_q,dvdx,dudy,νSmag_q,dLudy,dLvdx = preallocate_q_variables()

    RKa = Numtype.([1/6,1/3,1/3,1/6])
    RKb = Numtype.([.5,.5,1.])

    # feedback and output
    t0 = feedback_ini()
    ncs, iout = output_nc_ini(u,v,η)
    scripts_output()
    nans_detected = false

    t = 0           # model time
    for i = 1:nt

        u1[:],v1[:],η1[:] = u,v,η

        for rki = 1:4
            rhs!(du,dv,dη,u1,v1,η1,Fx,f_q,
                dudx,dvdy,dvdx,dudy,dpdx,dpdy,
                p,KEu,KEv,
                h,h_u,h_v,h_q,U,V,U_v,V_u,
                adv_u,adv_v,q,q_u,q_v,
                Lu1,Lu2,Lv1,Lv2,dLudx,dLudy,dLvdx,dLvdy,
                shear,νSmag,νSmag_q)


            if rki < 4
                u1[:] = u + RKb[rki]*Δt*du
                v1[:] = v + RKb[rki]*Δt*dv
                η1[:] = η + RKb[rki]*Δt*dη
            end

            u0 += RKa[rki]*Δt*du
            v0 += RKa[rki]*Δt*dv
            η0 += RKa[rki]*Δt*dη
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

function time_integration(u::AbstractMatrix,v::AbstractMatrix,η::AbstractMatrix)

    # FORCING
    Fx = channel_wind()
    f_u,f_v,f_q = beta_plane()
    H = seamount()

    # PREALLOCATE
    du,u0,u1,dpdx,U,V_u,h_u,q_u,adv_u,Lu1,Lu2 = preallocate_u_vars(u)
    dv,v0,v1,dpdy,V,U_v,h_v,q_v,adv_v,Lv1,Lv2 = preallocate_v_vars(v)
    dη,η0,η1,dudx,dvdy,dUdx,dVdy,h,KEu,KEv,p,νSmag,dLudx,dLvdy,shear = preallocate_T_variables(η)
    q,h_q,dvdx,dudy,νSmag_q,dLudy,dLvdx = preallocate_q_variables()

    # propagate initial conditions
    u0[:,:],v0[:,:],η0[:,:] = u,v,η

    RKa = Numtype.([1/6,1/3,1/3,1/6])
    RKb = Numtype.([.5,.5,1.])

    # feedback and output
    t0,progrtxt = feedback_ini()
    ncs, iout = output_ini(u,v,η)
    nans_detected = false

    t = 0           # model time
    for i = 1:nt

        u1[:,:],v1[:,:],η1[:,:] = u,v,η

        for rki = 1:4
            rhs!(du,dv,dη,u1,v1,η1,Fx,f_q,H,
                dudx,dvdy,dvdx,dudy,dpdx,dpdy,
                p,KEu,KEv,dUdx,dVdy,
                h,h_u,h_v,h_q,U,V,U_v,V_u,
                adv_u,adv_v,q,q_u,q_v,
                Lu1,Lu2,Lv1,Lv2,dLudx,dLudy,dLvdx,dLvdy,
                shear,νSmag,νSmag_q)


            if rki < 4
                @views u1[:,:] .= u .+ RKb[rki]*Δt*du
                @views v1[:,:] .= v .+ RKb[rki]*Δt*dv
                @views η1[:,:] .= η .+ RKb[rki]*Δt*dη
            end

            @views u0 .+= RKa[rki]*Δt*du
            @views v0 .+= RKa[rki]*Δt*dv
            @views η0 .+= RKa[rki]*Δt*dη
        end

        u[:,:],v[:,:],η[:,:] = u0,v0,η0
        t += dtint

        # feedback and output
        t0,nans_detected = feedback(u,v,η,i,t0,nt,nans_detected,progrtxt)
        ncs,iout = output_nc(ncs,u,v,η,i,iout)
    end

    # feeback and output
    feedback_end(progrtxt,t0)
    output_close(ncs,progrtxt)

    return u,v,η
end

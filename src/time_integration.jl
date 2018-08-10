function time_integration(u::AbstractMatrix,v::AbstractMatrix,η::AbstractMatrix)

    # FORCING
    Fx = wind()
    f_q = beta_plane()
    H = topography()

    # add halo with ghost point copy
    u,v,η = add_halo(u,v,η)

    # PREALLOCATE
    du,u0,u1,u²,KEu,dudx,dudy,Lu = preallocate_u_vars()
    dv,v0,v1,v²,KEv,dvdy,dvdx,Lv = preallocate_v_vars()
    dη,η0,η1,p,h,h_q,q,dpdx,h_u,U,dpdy,h_v,V,dUdx,dVdy,q_v,adv_v,U_v,q_u,adv_u,V_u = preallocate_T_variables()

    # propagate initial conditions
    u0 .= u
    v0 .= v
    η0 .= η

    RKa = Numtype.([1/6,1/3,1/3,1/6])
    RKb = Numtype.([.5,.5,1.])

    # feedback and output
    t0,progrtxt = feedback_ini()
    ncs, iout = output_ini(u,v,η)
    nans_detected = false

    t = 0           # model time
    for i = 1:nt

        ghost_points!(u,v,η)
        u1 .= u
        v1 .= v
        η1 .= η

        for rki = 1:4
            if rki > 1
                ghost_points!(u1,v1,η1)
            end

            rhs!(du,dv,dη,u1,v1,η1,Fx,f_q,H,
                dudx,dvdy,dvdx,dudy,dpdx,dpdy,
                p,u²,v²,KEu,KEv,dUdx,dVdy,
                h,h_u,h_v,h_q,U,V,U_v,V_u,
                adv_u,adv_v,q,q_u,q_v,
                Lu,Lv)
                #Lu1,Lu2,Lv1,Lv2,dLudx,dLudy,dLvdx,dLvdy,
                #shear,νSmag,νSmag_q)


            if rki < 4
                u1 .= u .+ RKb[rki]*Δt*du
                v1 .= v .+ RKb[rki]*Δt*dv
                η1 .= η .+ RKb[rki]*Δt*dη
            end

            u0 .+= RKa[rki]*Δt*du
            v0 .+= RKa[rki]*Δt*dv
            η0 .+= RKa[rki]*Δt*dη
        end

        u .= u0
        v .= v0
        η .= η0
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

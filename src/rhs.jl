function rhs(du,dv,dη,u,v,η,Fx,f_q,
             dpdx,dpdy,dLu,dLu2,dLv,dLv2,dudx,dvdy,
             p,KEu,KEv,
             h,h_u,h_v,U,V,U_v,V_u,
             q,dvdx,dudy,h_q,q_u,q_v,
             adv_u,adv_v)

    h[:] = η+H

    U[:] = u.*ITu(h_u,h)
    V[:] = v.*ITv(h_v,h)

    # Bernoulli potential
    p[:] = .5*(IuT(KEu,u.^2) + IvT(KEv,v.^2)) + g*h

    # Potential vorticity
    q[:] = (f_q + Gvx(dvdx,v) - Guy(dudy,u)) ./ ITq(h_q,h)

    # Sadourny, 1975 enstrophy conserving scheme
    adv_u[:] = Iqu(q_u,q) .* Ivu(V_u,V)
    adv_v[:] = -Iqv(q_v,q) .* Iuv(U_v,U)

    du[:] = adv_u - GTx(dpdx,p) + ν*Lu(dLu,u) + Fx
    dv[:] = adv_v - GTy(dpdy,p) + ν*Lv(dLv,v)
    dη[:] = -H*(Gux(dudx,u) + Gvy(dvdy,v))

    return du,dv,dη
end

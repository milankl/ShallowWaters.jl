function rhs(du,dv,dη,u,v,η,Fx,f_q,
             dpdx,dpdy,dLu,dLu2,dLv,dLv2,dudx,dvdy,
             p,KEu,KEv,
             h,h_u,h_v,U,V,U_v,V_u,
             q,dvdx,dudy,h_q,q_u,q_v,
             adv_u,adv_v,
             shear,νSmag,νSmag_u,νSmag_v)

    h[:] = η+H
    ITu(h_u,h)
    ITv(h_v,h)
    ITq(h_q,h)

    U[:] = u.*h_u
    V[:] = v.*h_v

    Gux(dudx,u)
    Gvy(dvdy,v)
    Gvx(dvdx,v)
    Guy(dudy,u)

    # Bernoulli potential
    IuT(KEu,u.^2)
    IvT(KEv,v.^2)
    p[:] = .5*(KEu + KEv) + g*h
    GTx(dpdx,p)
    GTy(dpdy,p)

    # Potential vorticity
    q[:] = (f_q + dvdx - dudy) ./ h_q

    # Sadourny, 1975 enstrophy conserving scheme
    Iqu(q_u,q)
    Iqv(q_v,q)
    Ivu(V_u,V)
    Iuv(U_v,U)
    adv_u[:] = q_u.*V_u
    adv_v[:] = -q_v.*U_v

    # Smagorinsky-like biharmonic diffusion
    #νSmag[:] = cSmag*sqrt.((dudx-dvdy).^2 + ITq(?,(dudy + dvdx).^2))
    #νSmag_u[:] = ITu(νSmag_u,νSmag)
    #νSmag_v[:] = ITv(νSmag_v,νSmag)
    Lu(dLu,u)
    Lv(dLv,v)

    # adding the terms
    du[:] = adv_u - dpdx + ν*dLu + Fx
    dv[:] = adv_v - dpdy + ν*dLv
    dη[:] = -H*(dudx + dvdy)

    return du,dv,dη
end

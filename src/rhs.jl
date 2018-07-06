function rhs!(du,dv,dη,u,v,η,Fx,f_q,
             dudx,dvdy,dvdx,dudy,dpdx,dpdy,
             p,KEu,KEv,
             h,h_u,h_v,h_q,U,V,U_v,V_u,
             adv_u,adv_v,q,q_u,q_v,
             Lu1,Lu2,Lv1,Lv2,dLudx,dLudy,dLvdx,dLvdy,
             shear,νSmag,νSmag_q)

    h[:] = η+H
    ITu!(h_u,h)
    ITv!(h_v,h)
    ITq!(h_q,h)

    U[:] = u.*h_u
    V[:] = v.*h_v

    Gux!(dudx,u)
    Gvy!(dvdy,v)
    Gvx!(dvdx,v)
    Guy!(dudy,u)

    # Bernoulli potential
    IuT!(KEu,u.^2)
    IvT!(KEv,v.^2)
    p[:] = one_half*(KEu + KEv) + g*η
    GTx!(dpdx,p)
    GTy!(dpdy,p)

    # Potential vorticity
    q[:] = (f_q + dvdx - dudy) ./ h_q

    # Sadourny, 1975 enstrophy conserving scheme
    Iqu!(q_u,q)
    Iqv!(q_v,q)
    Ivu!(V_u,V)
    Iuv!(U_v,U)
    adv_u[:] = q_u.*V_u
    adv_v[:] = -q_v.*U_v

    #= Smagorinsky-like biharmonic diffusion
    Lu1 + Lu2 = dx[ νSmag dx(L(u))] + dy[ νSmag dy(L(u))]
    Lv1 + Lv2 = dx[ νSmag dx(L(v))] + dy[ νSmag dy(L(v))]
    =#
    IqT!(shear,(dudy + dvdx).^2)
    νSmag[:] = cSmag*sqrt.((dudx-dvdy).^2 + shear)
    ITq!(νSmag_q,νSmag)
    ∇²u!(Lu1,u)
    ∇²v!(Lv1,v)
    Gux!(dLudx,Lu1)
    Guy!(dLudy,Lu1)
    Gvx!(dLvdx,Lv1)
    Gvy!(dLvdy,Lv1)
    GTx!(Lu1,νSmag.*dLudx)        # reuse intermediate variable Lu1
    Gqy!(Lu2,νSmag_q.*dLudy)
    GTy!(Lv1,νSmag.*dLvdy)
    Gqx!(Lv2,νSmag_q.*dLvdx)

    # adding the terms
    du[:] = adv_u - dpdx + Lu1+Lu2 + Fx
    dv[:] = adv_v - dpdy + Lv1+Lv2
    dη[:] = -H*(dudx + dvdy)
end

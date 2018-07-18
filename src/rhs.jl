function rhs!(du,dv,dη,u,v,η,Fx,f_q,H,
             dudx,dvdy,dvdx,dudy,dpdx,dpdy,
             p,KEu,KEv,dUdx,dVdy,
             h,h_u,h_v,h_q,U,V,U_v,V_u,
             adv_u,adv_v,q,q_u,q_v,
             Lu1,Lu2,Lv1,Lv2,dLudx,dLudy,dLvdx,dLvdy,
             shear,νSmag,νSmag_q)

    @views h[:,:] .= η .+ H
    ITu!(h_u,h)     # use η instead of h preliminary
    ITv!(h_v,h)
    ITq!(h_q,h)

    @views U[:,:] .= u.*h_u
    @views V[:,:] .= v.*h_v

    Gux!(dudx,u)
    Gvy!(dvdy,v)
    Gvx!(dvdx,v)
    Guy!(dudy,u)

    Gux!(dUdx,U)
    Gvy!(dVdy,V)

    # Bernoulli potential
    @views Lu1[:,:] .= u.^2         # reuse Lu1, Lv1
    @views Lv1[:,:] .= v.^2
    IuT!(KEu,Lu1)
    IvT!(KEv,Lv1)
    @views p[:,:] .= one_half*(KEu .+ KEv) .+ g*η
    GTx!(dpdx,p)
    GTy!(dpdy,p)

    # Potential vorticity
    @views q[:,:] .= (f_q .+ dvdx .- dudy) ./ h_q

    # Sadourny, 1975 enstrophy conserving scheme
    Iqu!(q_u,q)
    Iqv!(q_v,q)
    Ivu!(V_u,V)
    Iuv!(U_v,U)
    @views adv_u[:,:] .= q_u.*V_u
    @views adv_v[:,:] .= -q_v.*U_v

    #= Smagorinsky-like biharmonic diffusion
    Lu1 + Lu2 = dx[ νSmag dx(L(u))] + dy[ νSmag dy(L(u))]
    Lv1 + Lv2 = dx[ νSmag dx(L(v))] + dy[ νSmag dy(L(v))]
    =#
    @views νSmag_q[:,:] .= (dudy .+ dvdx).^2     # reuse variable νSmag_q
    IqT!(shear,νSmag_q)
    @views νSmag[:,:] .= cSmag*sqrt.((dudx .- dvdy).^2 .+ shear)
    ITq!(νSmag_q,νSmag)
    ∇²u!(Lu1,u)
    ∇²v!(Lv1,v)
    Gux!(dLudx,Lu1)
    Guy!(dLudy,Lu1)
    Gvx!(dLvdx,Lv1)
    Gvy!(dLvdy,Lv1)
    @views dLudx[:,:] .= νSmag.*dLudx
    @views dLudy[:,:] .= νSmag_q.*dLudy
    @views dLvdy[:,:] .= νSmag.*dLvdy
    @views dLvdx[:,:] .= νSmag_q.*dLvdx
    GTx!(Lu1,dLudx)        # reuse intermediate variable Lu1
    Gqy!(Lu2,dLudy)
    GTy!(Lv1,dLvdy)
    Gqx!(Lv2,dLvdx)

    # adding the terms
    @views du[:,:] .= adv_u .- dpdx .+ Lu1.+Lu2 .+ Fx
    @views dv[:,:] .= adv_v .- dpdy .+ Lv1.+Lv2
    @views dη[:,:] .= -(dUdx .+ dVdy)
end

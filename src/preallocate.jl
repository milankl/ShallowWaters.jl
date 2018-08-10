function preallocate_u_vars()

    # with full halo
    du = zeros(Numtype,nux+2*halo,nuy+2*halo)
    u0 = zeros(du)
    u1 = zeros(du)
    u² = zeros(du)

    # one less in x-direction
    KEu = zeros(Numtype,nux+2*halo-1,nuy+2*halo)
    dudx = zeros(KEu)

    # one less in y-direction
    dudy = zeros(Numtype,nux+2*halo,nuy+2*halo-1)

    # two less in x-direction one less in y-direction
    U_v = zeros(Numtype,nux+2*halo-2,nuy+2*halo-1)

    # two less in both directions
    Lu = zeros(Numtype,nux+2*halo-2,nuy+2*halo-2)


    #TODO
    #Lu1 = zeros(u)
    #Lu2 = zeros(u)

    return du,u0,u1,u²,KEu,dudx,dudy,Lu
end

function preallocate_v_vars()

    # with full halo
    dv = zeros(Numtype,nvx+2*halo,nvy+2*halo)
    v0 = zeros(dv)
    v1 = zeros(dv)
    v² = zeros(dv)

    # one less in y-direction
    KEv = zeros(Numtype,nvx+2*halo,nvy+2*halo-1)
    dvdy = zeros(KEv)

    # one less in x-direction
    dvdx = zeros(Numtype,nvx+2*halo-1,nvy+2*halo)

    # two less in both directions
    Lv = zeros(Numtype,nvx+2*halo-2,nvy+2*halo-2)

    #TODO
    #Lv1 = zeros(v)
    #Lv2 = zeros(v)

    return dv,v0,v1,v²,KEv,dvdy,dvdx,Lv
end

function preallocate_T_variables()

    # full halo
    dη = zeros(Numtype,nx+2,ny+2)  # halo for η is always 1
    η0 = zeros(dη)
    η1 = zeros(dη)
    p = zeros(dη)
    h = zeros(dη)

    # one less in x and y-direction
    h_q = zeros(Numtype,nx+1,ny+1)
    q = zeros(h_q)

    # one less in x direction
    dpdx = zeros(Numtype,nx+1,ny+2)
    h_u = zeros(dpdx)
    U = zeros(dpdx)

    # one less in y-direction
    dpdy = zeros(Numtype,nx+2,ny+1)
    h_v = zeros(dpdy)
    V = zeros(dpdy)

    # two less in x-direction
    dUdx = zeros(Numtype,nx,ny+2)

    # two less in y-direction
    dVdy = zeros(Numtype,nx+2,ny)

    # two less in x direction, one less in y
    q_v = zeros(Numtype,nx,ny+1)
    adv_v = zeros(q_v)
    U_v = zeros(q_v)

    # two less in y direction, one less in x
    q_u = zeros(Numtype,nx+1,ny)
    adv_u = zeros(q_u)
    V_u = zeros(q_u)

    #TODO
    #νSmag = zeros(η)
    #dLudx = zeros(η)
    #dLvdy = zeros(η)
    #shear = zeros(η)

    return dη,η0,η1,p,h,h_q,q,dpdx,h_u,U,dpdy,h_v,V,dUdx,dVdy,q_v,adv_v,U_v,q_u,adv_u,V_u
end

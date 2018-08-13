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

    return du,u0,u1,u²,KEu,dudx,dudy
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

    return dv,v0,v1,v²,KEv,dvdy,dvdx
end

function preallocate_T_variables()

    # full halo
    dη = zeros(Numtype,nx+2,ny+2)  # halo for η is always 1
    η0 = zeros(dη)
    η1 = zeros(dη)
    p = zeros(dη)
    h = zeros(dη)

    # one less in both directions
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
    qhu = zeros(q_v)
    U_v = zeros(q_v)

    # two less in y direction, one less in x
    q_u = zeros(Numtype,nx+1,ny)
    qhv = zeros(q_u)
    V_u = zeros(q_u)

    return dη,η0,η1,p,h,h_q,q,dpdx,h_u,U,dpdy,h_v,V,dUdx,dVdy,q_v,qhu,U_v,q_u,qhv,V_u
end

function preallocate_Sadourny()
    #TODO
end

function preallocate_Smagorinsky()
    # on the T-grid including halo
    DT = zeros(Numtype,nx+2,ny+2)
    DS = zeros(DT)
    νSmag = zeros(DT)

    # DS_q has same size as dvdx
    DS_q = zeros(Numtype,nvx+2*halo-1,nvy+2*halo)

    # one less in both directions, the q-grid
    νSmag_q = zeros(Numtype,nx+1,ny+1)
    S12 = zeros(νSmag_q)
    S21 = zeros(νSmag_q)

    # Laplace operator: two less in both directions
    Lu = zeros(Numtype,nux+2*halo-2,nuy+2*halo-2)
    Lv = zeros(Numtype,nvx+2*halo-2,nvy+2*halo-2)

    # 3 less in x, two less in y
    dLudx = zeros(Numtype,nux+2*halo-3,nuy+2*halo-2)
    S11 = zeros(dLudx)

    # 4 less in x, two less in y
    LLu1 = zeros(Numtype,nux+2*halo-4,nuy+2*halo-2)

    # 2 less in x, 3 less in y
    dLudy = zeros(Numtype,nux+2*halo-2,nuy+2*halo-3)

    LLu2 = zeros(Numtype,nx+1,ny)

    # three less in x, two less in y
    dLvdx = zeros(Numtype,nvx+2*halo-3,nvy+2*halo-2)

    # 4 less in x, two less in y
    LLv1 = zeros(Numtype,nx,ny+1)

    # two less in x, three less in y
    dLvdy = zeros(Numtype,nvx+2*halo-2,nvy+2*halo-3)
    S22 = zeros(dLvdy)

    # 2 less in x, 4 less in y
    LLv2 = zeros(Numtype,nvx+2*halo-2,nvy+2*halo-4)

    return DT,DS,DS_q,νSmag,νSmag_q,Lu,Lv,dLudx,dLudy,dLvdx,dLvdy,S11,S12,S21,S22,LLu1,LLu2,LLv1,LLv2
end

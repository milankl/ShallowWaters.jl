function preallocate_u_vars()

    # with full halo
    du = zeros(Numtype,nux+2*halo,nuy+2*halo)
    u0 = zeros(du)
    u1 = zeros(du)

    # derivative: one less in x or y direction
    dudx = zeros(Numtype,nux+2*halo-1,nuy+2*halo)
    dudy = zeros(Numtype,nux+2*halo,nuy+2*halo-1)

    return du,u0,u1,dudx,dudy
end

function preallocate_v_vars()

    # with full halo
    dv = zeros(Numtype,nvx+2*halo,nvy+2*halo)
    v0 = zeros(dv)
    v1 = zeros(dv)

    # derivative: one less in x or y direction
    dvdx = zeros(Numtype,nvx+2*halo-1,nvy+2*halo)
    dvdy = zeros(Numtype,nvx+2*halo,nvy+2*halo-1)

    return dv,v0,v1,dvdx,dvdy
end

function preallocate_η_vars()

    dη = zeros(Numtype,nx+2*haloη,ny+2*haloη)
    η0 = zeros(dη)
    η1 = zeros(dη)
    h = zeros(dη)

    return dη,η0,η1,h
end

function preallocate_continuity()

    # interpolation: one less in x-direction
    h_u = zeros(Numtype,nx+2*haloη-1,ny+2*haloη)
    U = zeros(h_u)

    # interpolation: one less in y-direction
    h_v = zeros(Numtype,nx+2*haloη,ny+2*haloη-1)
    V = zeros(h_v)

    # Derivatives: two less in x- or y-direction
    dUdx = zeros(Numtype,nx+2*haloη-2,ny+2*haloη)
    dVdy = zeros(Numtype,nx+2*haloη,ny+2*haloη-2)

    return h_u,U,h_v,V,dUdx,dVdy
end

function preallocate_Sadourny()

    # interpolation from h: ones less in both directions
    h_q = zeros(Numtype,nx+2*haloη-1,ny+2*haloη-1)
    q = zeros(h_q)

    # two less in x direction, one less in y
    q_v = zeros(Numtype,nx+2*haloη-2,ny+2*haloη-1)
    qhu = zeros(q_v)
    U_v = zeros(q_v)

    # two less in y direction, one less in x
    q_u = zeros(Numtype,nx+2*haloη-1,ny+2*haloη-2)
    qhv = zeros(q_u)
    V_u = zeros(q_u)

    return h_q,q,q_v,qhu,U_v,q_u,qhv,V_u
end

function preallocate_Bernoulli()

    u² = zeros(Numtype,nux+2*halo,nuy+2*halo)
    v² = zeros(Numtype,nvx+2*halo,nvy+2*halo)

    KEu = zeros(Numtype,nux+2*halo-1,nuy+2*halo)
    KEv = zeros(Numtype,nvx+2*halo,nvy+2*halo-1)

    p = zeros(Numtype,nx+2*haloη,ny+2*haloη)
    dpdx = zeros(Numtype,nx+2*haloη-1,ny+2*haloη)
    dpdy = zeros(Numtype,nx+2*haloη,ny+2*haloη-1)

    return u²,v²,KEu,KEv,p,dpdx,dpdy
end

function preallocate_Smagorinsky()
    # on the T-grid including halo
    DT = zeros(Numtype,nx+2*haloη,ny+2*haloη)
    DS = zeros(DT)
    νSmag = zeros(DT)

    # DS_q has same size as dvdx
    DS_q = zeros(Numtype,nvx+2*halo-1,nvy+2*halo)

    # one less in both directions, the q-grid
    νSmag_q = zeros(Numtype,nx+2*haloη-1,ny+2*haloη-1)
    S12 = zeros(νSmag_q)
    S21 = zeros(νSmag_q)

    S11 = zeros(Numtype,nux+2*halo-3,nuy+2*halo-2)
    S22 = zeros(Numtype,nvx+2*halo-2,nvy+2*halo-3)

    LLu1 = zeros(Numtype,nux+2*halo-4,nuy+2*halo-2)
    LLu2 = zeros(Numtype,nx+1,ny)

    LLv1 = zeros(Numtype,nx,ny+1)
    LLv2 = zeros(Numtype,nvx+2*halo-2,nvy+2*halo-4)

    return DT,DS,DS_q,νSmag,νSmag_q,S11,S12,S21,S22,LLu1,LLu2,LLv1,LLv2
end

function preallocate_Laplace()
    # two less in both directions
    Lu = zeros(Numtype,nux+2*halo-2,nuy+2*halo-2)
    Lv = zeros(Numtype,nvx+2*halo-2,nvy+2*halo-2)

    # Derivatives of Lu,Lv
    dLudx = zeros(Numtype,nux+2*halo-3,nuy+2*halo-2)
    dLudy = zeros(Numtype,nux+2*halo-2,nuy+2*halo-3)
    dLvdx = zeros(Numtype,nvx+2*halo-3,nvy+2*halo-2)
    dLvdy = zeros(Numtype,nvx+2*halo-2,nvy+2*halo-3)

    return Lu,Lv,dLudx,dLudy,dLvdx,dLvdy
end

function preallocate_Laplace()
    # two less in both directions
    Lu = zeros(Numtype,nux+2*halo-2,nuy+2*halo-2)
    Lv = zeros(Numtype,nvx+2*halo-2,nvy+2*halo-2)

    # Derivatives of Lu,Lv
    dLudx = zeros(Numtype,nux+2*halo-3,nuy+2*halo-2)
    dLudy = zeros(Numtype,nux+2*halo-2,nuy+2*halo-3)
    dLvdx = zeros(Numtype,nvx+2*halo-3,nvy+2*halo-2)
    dLvdy = zeros(Numtype,nvx+2*halo-2,nvy+2*halo-3)

    return Lu,Lv,dLudx,dLudy,dLvdx,dLvdy
end

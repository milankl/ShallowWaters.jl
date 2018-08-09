function preallocate_u_vars()

    # with full halo
    du = zeros(Numtype,nux+2*halo,nuy+2*halo)
    u0 = zeros(du)
    u1 = zeros(du)
    u² = zeros(du)

    # one less in x-direction
    U = zeros(Numtype,nux+2*halo-1,nuy+2*halo)
    h_u = zeros(U)
    KEu = zeros(U)
    dudx = zeros(U)

    # one less in y-direction
    dudy = zeros(Numtype,nux+2*halo,nuy+2*halo-1)

    # two less in x-direction
    dUdx = zeros(Numtype,nux+2*halo-2,nuy+2*halo)


    #TODO
    V_u = zeros(u)

    q_u = zeros(u)
    adv_u = zeros(u)

    Lu1 = zeros(u)
    Lu2 = zeros(u)

    return
end

function preallocate_v_vars()

    # with full halo
    dv = zeros(Numtype,nvx+2*halo,nvy+2*halo)
    v0 = zeros(dv)
    v1 = zeros(dv)
    v² = zeros(dv)

    # one less in y-direction
    V = zeros(Numtype,nvx+2*halo,nvy+2*halo-1)
    h_v = zeros(V)
    KEv = zeros(V)
    dvdy = zeros(V)

    # one less in x-direction
    dvdx = zeros(Numtype,nvx+2*halo-1,nvy+2*halo)

    # two less in y-direction
    dVdy = zeros(Numtype,nvx+2*halo,nvy+2*halo-2)

    #TODO
    dpdy = zeros(v)

    U_v = zeros(v)

    q_v = zeros(v)
    adv_v = zeros(v)

    Lv1 = zeros(v)
    Lv2 = zeros(v)

    return
end

function preallocate_T_variables()

    # full halo
    dη = zeros(Numtype,nx+2,nvy+2)  # halo for η is always 1
    η0 = zeros(dη)
    η1 = zeros(dη)
    p = zeros(dη)
    h = zeros(dη)

    # one less in x and y-direction
    h_q = zeros(Numtype,nx+1,ny+1)

    # one less in x direction
    dpdx = zeros(Numtype,nx+1,ny+2)

    # one less in y-direction
    dpdy = zeros(Numtype,nx+2,ny+1)

    #TODO
    νSmag = zeros(η)
    dLudx = zeros(η)
    dLvdy = zeros(η)
    shear = zeros(η)

    return
end

function preallocate_u_vars(u::AbstractMatrix)

    # initialise with zeros
    du = zeros(u)
    u0 = zeros(u)
    u1 = zeros(u)

    dpdx = zeros(u)
    U = zeros(u)
    V_u = zeros(u)
    h_u = zeros(u)
    q_u = zeros(u)
    adv_u = zeros(u)

    Lu1 = zeros(u)
    Lu2 = zeros(u)

    return du,u0,u1,dpdx,U,V_u,h_u,q_u,adv_u,Lu1,Lu2
end

function preallocate_v_vars(v::AbstractMatrix)

    # initialise with zeros
    dv = zeros(v)
    v0 = zeros(v)
    v1 = zeros(v)

    dpdy = zeros(v)
    V = zeros(v)
    U_v = zeros(v)
    h_v = zeros(v)
    q_v = zeros(v)
    adv_v = zeros(v)

    Lv1 = zeros(v)
    Lv2 = zeros(v)

    return dv,v0,v1,dpdy,V,U_v,h_v,q_v,adv_v,Lv1,Lv2
end

function preallocate_T_variables(η::AbstractMatrix)

    # initialise with zeros
    dη = zeros(η)
    η0 = zeros(η)
    η1 = zeros(η)

    dudx = zeros(η)
    dvdy = zeros(η)
    dUdx = zeros(η)
    dVdy = zeros(η)
    h = zeros(η)
    KEu = zeros(η)
    KEv = zeros(η)
    p = zeros(η)

    νSmag = zeros(η)
    dLudx = zeros(η)
    dLvdy = zeros(η)
    shear = zeros(η)

    return dη,η0,η1,dudx,dvdy,dUdx,dVdy,h,KEu,KEv,p,νSmag,dLudx,dLvdy,shear
end

function preallocate_q_variables()

    # initialise with zeros
    q = zeros(Numtype,nqx,nqy)
    h_q = zeros(Numtype,nqx,nqy)
    dvdx = zeros(Numtype,nqx,nqy)
    dudy = zeros(Numtype,nqx,nqy)

    νSmag_q = zeros(Numtype,nqx,nqy)
    dLudy = zeros(Numtype,nqx,nqy)
    dLvdx = zeros(Numtype,nqx,nqy)

    return q,h_q,dvdx,dudy,νSmag_q,dLudy,dLvdx
end

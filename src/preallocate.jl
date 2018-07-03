function preallocate_u_vars(u)

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

    dLu = zeros(u)

    return du,u0,u1,dpdx,U,V_u,h_u,q_u,adv_u,dLu
end

function preallocate_v_vars(v)

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

    dLv = zeros(v)

    return dv,v0,v1,dpdy,V,U_v,h_v,q_v,adv_v,dLv
end

function preallocate_T_variables(η)

    # initialise with zeros
    dη = zeros(η)
    η0 = zeros(η)
    η1 = zeros(η)

    dudx = zeros(η)
    dvdy = zeros(η)
    h = zeros(η)
    KEu = zeros(η)
    KEv = zeros(η)
    p = zeros(η)

    return dη,η0,η1,dudx,dvdy,h,KEu,KEv,p
end

function preallocate_q_variables()

    # initialise with zeros
    q = zeros(nqx,nqy)
    h_q = zeros(nqx,nqy)
    dvdx = zeros(nqx,nqy)
    dudy = zeros(nqx,nqy)

    return q,h_q,dvdx,dudy
end

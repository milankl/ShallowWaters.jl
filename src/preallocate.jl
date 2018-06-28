function preallocate_variables(u,v,η)

    # initialise with zeros
    v_u = zeros(u)
    u_v = zeros(v)

    dηdx = zeros(u)
    dηdy = zeros(v)

    dLu = zeros(u)
    dLv = zeros(v)

    dudx = zeros(η)
    dvdy = zeros(η)

    return v_u,u_v,dηdx,dηdy,dLu,dLv,dudx,dvdy
end

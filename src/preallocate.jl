function preallocate_u_vars(u)

    # initialise with zeros
    v_u = zeros(u)
    dηdx = zeros(u)
    dLu = zeros(u)

    return v_u,dηdx,dLu
end

function preallocate_v_vars(v)

    # initialise with zeros
    u_v = zeros(v)
    dηdy = zeros(v)
    dLv = zeros(v)

    return u_v,dηdy,dLv
end

function preallocate_T_variables(η)

    # initialise with zeros
    dudx = zeros(η)
    dvdy = zeros(η)

    return dudx,dvdy
end

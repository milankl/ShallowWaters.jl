function preallocate_u_vars(u)

    # initialise with zeros
    v_u = zeros(u)
    dηdx = zeros(u)
    dLu = zeros(u)
    u0 = zeros(u)
    u1 = zeros(u)
    du = zeros(u)

    return du,u0,u1,v_u,dηdx,dLu
end

function preallocate_v_vars(v)

    # initialise with zeros
    u_v = zeros(v)
    dηdy = zeros(v)
    dLv = zeros(v)
    v0 = zeros(v)
    v1 = zeros(v)
    dv = zeros(v)

    return dv,v0,v1,u_v,dηdy,dLv
end

function preallocate_T_variables(η)

    # initialise with zeros
    dudx = zeros(η)
    dvdy = zeros(η)
    η0 = zeros(η)
    η1 = zeros(η)
    dη = zeros(η)

    return dη,η0,η1,dudx,dvdy
end

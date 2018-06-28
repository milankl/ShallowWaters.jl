function rhs(u,v,η,v_u,u_v,dηdx,dηdy,dLu,dLv,dudx,dvdy)

    du[:] = f_0*Ivu(v_u,v) - g*GTx(dηdx,η) + ν*Lu(dLu,u) + Fx
    dv[:] = -f_0*Iuv(u_v,u) - g*GTy(dηdy,η) + ν*Lv(dLv,v)
    dη[:] = -H*(Gux(dudx,u) + Gvy(dvdy,v))

    return du,dv,dη
end

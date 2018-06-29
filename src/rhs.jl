function rhs(du,dv,dη,u,v,η,Fx,f_u,f_v,v_u,u_v,dηdx,dηdy,dLu,dLv,dudx,dvdy)

    du[:] = f_u.*Ivu(v_u,v) - gg*GTx(dηdx,η) + ν*Lu(dLu,u) + Fx
    dv[:] = -f_v.*Iuv(u_v,u) - gg*GTy(dηdy,η) + ν*Lv(dLv,v)
    dη[:] = -HH*(Gux(dudx,u) + Gvy(dvdy,v))

    return du,dv,dη
end

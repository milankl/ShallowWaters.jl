function initial_conditions()
    u_0 = zeros(Nux,Nuy)
    v_0 = zeros(Nvx,Nvy)
    eta_0 = zeros(NTx,NTy)
    return u_0,v_0,eta_0
end

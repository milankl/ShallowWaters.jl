function initial_conditions()
    u = zeros(Numtype,nux,nuy)
    v = zeros(Numtype,nvx,nvy)
    η = zeros(Numtype,nx,ny)
    return u,v,η
end

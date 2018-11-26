"""Returns preallocated variables of different size that derive from u."""
function preallocate_u_vars()
    # with full halo
    du = zeros(Numtype,nux+2*halo,nuy+2*halo)
    u0 = zero(du)
    u1 = zero(du)

    # derivative: one less in x or y direction
    dudx = zeros(Numtype,nux+2*halo-1,nuy+2*halo)
    dudy = zeros(Numtype,nux+2*halo,nuy+2*halo-1)

    return du,u0,u1,dudx,dudy
end

"""Returns preallocated variables of different size that derive from v."""
function preallocate_v_vars()
    # with full halo
    dv = zeros(Numtype,nvx+2*halo,nvy+2*halo)
    v0 = zero(dv)
    v1 = zero(dv)

    # derivative: one less in x or y direction
    dvdx = zeros(Numtype,nvx+2*halo-1,nvy+2*halo)
    dvdy = zeros(Numtype,nvx+2*halo,nvy+2*halo-1)

    return dv,v0,v1,dvdx,dvdy
end

"""Returns preallocated variables of different size that derive from η."""
function preallocate_η_vars()
    # with full halo
    dη = zeros(Numtype,nx+2*haloη,ny+2*haloη)
    η0 = zero(dη)
    η1 = zero(dη)
    h = zero(dη)

    return dη,η0,η1,h
end

"""Returns preallocated variables of different size that will be used in the contuinity equation."""
function preallocate_continuity()
    # interpolation: one less in x-direction
    h_u = zeros(Numtype,nx+2*haloη-1,ny+2*haloη)
    U = zero(h_u)

    # interpolation: one less in y-direction
    h_v = zeros(Numtype,nx+2*haloη,ny+2*haloη-1)
    V = zero(h_v)

    # Derivatives: two less in x- or y-direction
    dUdx = zeros(Numtype,nx+2*haloη-2,ny+2*haloη)
    dVdy = zeros(Numtype,nx+2*haloη,ny+2*haloη-2)

    return h_u,U,h_v,V,dUdx,dVdy
end

"""Returns preallocated variables of different size that will be used for PV and the Sadourny adv scheme."""
function preallocate_Sadourny()

    # interpolation from h: ones less in both directions
    h_q = zeros(Numtype,nx+2*haloη-1,ny+2*haloη-1)
    q = zero(h_q)

    # two less in x direction, one less in y
    q_v = zeros(Numtype,nx+2*haloη-2,ny+2*haloη-1)
    U_v = zero(q_v)

    # two less in y direction, one less in x
    q_u = zeros(Numtype,nx+2*haloη-1,ny+2*haloη-2)
    V_u = zero(q_u)

    #
    qhu = zeros(Numtype,nvx,nvy)
    qhv = zeros(Numtype,nux,nuy)

    return h_q,q,q_v,qhu,U_v,q_u,qhv,V_u
end

"""Returns preallocated variables of different size for the PV linear combinations in the Arakawa&Hsu advection scheme."""
function preallocate_ArakawaHsu()
    qα = zeros(Numtype,nx+2*haloη-2,ny+2*haloη-2)
    qβ = zeros(Numtype,nx+2*haloη-1,ny+2*haloη-2)
    qγ = zeros(Numtype,nx+2*haloη-1,ny+2*haloη-2)
    qδ = zeros(Numtype,nx+2*haloη-2,ny+2*haloη-2)
    return qα,qβ,qγ,qδ
end

"""Returns preallocated variables of different size that will be used for the Bernoulli potential."""
function preallocate_Bernoulli()
    # with full halo
    u² = zeros(Numtype,nux+2*halo,nuy+2*halo)
    v² = zeros(Numtype,nvx+2*halo,nvy+2*halo)

    KEu = zeros(Numtype,nux+2*halo-1,nuy+2*halo)
    KEv = zeros(Numtype,nvx+2*halo,nvy+2*halo-1)

    p = zeros(Numtype,nx+2*haloη,ny+2*haloη)
    dpdx = zeros(Numtype,nx+2*haloη-1,ny+2*haloη)
    dpdy = zeros(Numtype,nx+2*haloη,ny+2*haloη-1)

    return u²,v²,KEu,KEv,p,dpdx,dpdy
end

"""Returns preallocated variables of different size that will be used for the bottom drag term."""
function preallocate_bottomdrag()
    sqrtKE = zeros(Numtype,nx+2*haloη,ny+2*haloη)
    sqrtKE_u = zeros(Numtype,nx+2*haloη-1,ny+2*haloη)
    sqrtKE_v = zeros(Numtype,nx+2*haloη,ny+2*haloη-1)

    Bu = zeros(Numtype,nx+2*haloη-1,ny+2*haloη)
    Bv = zeros(Numtype,nx+2*haloη,ny+2*haloη-1)

    return sqrtKE,sqrtKE_u,sqrtKE_v,Bu,Bv
end

"""Returns preallocated variables of different size that will be used for the momentum diffusion term."""
function preallocate_Laplace()
    # two less in both directions
    Lu = zeros(Numtype,nux+2*halo-2,nuy+2*halo-2)
    Lv = zeros(Numtype,nvx+2*halo-2,nvy+2*halo-2)

    # Derivatives of Lu,Lv
    dLudx = zeros(Numtype,nux+2*halo-3,nuy+2*halo-2)
    dLudy = zeros(Numtype,nux+2*halo-2,nuy+2*halo-3)
    dLvdx = zeros(Numtype,nvx+2*halo-3,nvy+2*halo-2)
    dLvdy = zeros(Numtype,nvx+2*halo-2,nvy+2*halo-3)

    return Lu,Lv,dLudx,dLudy,dLvdx,dLvdy
end

"""Returns preallocated variables of different size that will be used in the contuinity equation."""
function preallocate_Smagorinsky()
    # on the η-grid including halo
    DT = zeros(Numtype,nx+2*haloη,ny+2*haloη)
    DS = zero(DT)
    νSmag = zero(DT)

    # DS_q has same size as dvdx
    DS_q = zeros(Numtype,nvx+2*halo-1,nvy+2*halo)

    # one less in both directions, the q-grid
    νSmag_q = zeros(Numtype,nx+2*haloη-1,ny+2*haloη-1)
    S12 = zero(νSmag_q)
    S21 = zero(νSmag_q)

    S11 = zeros(Numtype,nux+2*halo-3,nuy+2*halo-2)
    S22 = zeros(Numtype,nvx+2*halo-2,nvy+2*halo-3)

    LLu1 = zeros(Numtype,nux+2*halo-4,nuy+2*halo-2)
    LLu2 = zeros(Numtype,nx+1,ny)

    LLv1 = zeros(Numtype,nx,ny+1)
    LLv2 = zeros(Numtype,nvx+2*halo-2,nvy+2*halo-4)

    return DT,DS,DS_q,νSmag,νSmag_q,S11,S12,S21,S22,LLu1,LLu2,LLv1,LLv2
end

"""Returns preallocated variables of different size that will be used for the semi-Lagrangian tracer advection."""
function preallocate_semiLagrange()
    xd = zeros(Numtype,nx,ny)
    yd = zero(xd)

    um = zeros(Numtype,nux+2*halo,nuy+2*halo)
    vm = zeros(Numtype,nvx+2*halo,nvy+2*halo)

    # u on T-grid one less in x-direction
    u_T = zeros(Numtype,nux+2*halo-1,nuy+2*halo)
    um_T = zero(u_T)

    # v on T-grid one less in y-direction
    v_T = zeros(Numtype,nvx+2*halo,nvy+2*halo-1)
    vm_T = zero(v_T)

    uinterp = zeros(Numtype,nx,ny)
    vinterp = zero(uinterp)

    ssti = zeros(Numtype,nx+2*halosstx,ny+2*halossty)

    return xd,yd,um,vm,u_T,um_T,v_T,vm_T,uinterp,vinterp,ssti
end

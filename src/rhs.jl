"""Transit function to call either the linear or the nonlinear rhs."""
function rhs!(  u::AbstractMatrix,
                v::AbstractMatrix,
                η::AbstractMatrix,
                P::Parameter,
                C::Constants,
                G::Grid,
                Diag::DiagnosticVars,
                Forc::Forcing)

    if P.dynamics == "linear"
        rhs_linear!(u,v,η,P,C,G,Diag,Forc)
    else
        rhs_nonlinear!(u,v,η,P,C,G,Diag,Forc)
    end
end

"""Tendencies du,dv,dη of

        ∂u/∂t = qhv - ∂(1/2*(u²+v²) + gη)/∂x + Fx
        ∂v/∂t = -qhu - ∂(1/2*(u²+v²) + gη)/∂y + Fy
        ∂η/∂t =  -∂(uh)/∂x - ∂(vh)/∂y + γ(η_ref-η) + Fηt*Fη

the nonlinear shallow water equations."""
function rhs_nonlinear!(u::AbstractMatrix,
                        v::AbstractMatrix,
                        η::AbstractMatrix,
                        P::Parameter,
                        C::Constants,
                        G::Grid,
                        Diag::DiagnosticVars,
                        Forc::Forcing)

    @unpack h,h_u,h_v,U,V,dUdx,dVdy = Diag.VolumeFluxes
    @unpack H = Forc
    @unpack g = C
    @unpack ep = G

    # layer thickness
    thickness!(h,η,H)
    Ix!(h_u,h)
    Iy!(h_v,h)

    # mass or volume flux U,V = uh,vh
    Uflux!(U,u,h_u,ep)
    Vflux!(V,v,h_v)

    # divergence of mass flux
    ∂x!(dUdx,U)
    ∂y!(dVdy,V)

    if G.nstep_advcor == 0    # evaluate every RK substep
        advection_coriolis!(u,v,η,P,G,Diag,Forc)
    end

    # Bernoulli potential - recalculate for new η, KEu,KEv are only updated outside
    @unpack p,KEu,KEv,dpdx,dpdy = Diag.Bernoulli
    bernoulli!(p,KEu,KEv,η,g,ep)
    ∂x!(dpdx,p)
    ∂y!(dpdy,p)

    # Potential vorticity and advection thereof
    PVadvection!(P,G,Diag)

    # adding the terms
    @unpack du,dv,dη = Diag.Tendencies
    @unpack qhv,qhu = Diag.Vorticity
    @unpack Fx,Fy = Forc

    momentum_u!(du,qhv,dpdx,Fx,ep)
    momentum_v!(dv,qhu,dpdy,Fy)
    continuity!(dη,dUdx,dVdy)
end

"""Tendencies du,dv,dη of

        ∂u/∂t = gv - g∂η/∂x + Fx
        ∂v/∂t = -fu - g∂η/∂y
        ∂η/∂t =  -∂(uH)/∂x - ∂(vH)/∂y + γ(η_ref-η) + Fηt*Fη,

the linear shallow water equations."""
function rhs_linear!(   u::AbstractMatrix,
                        v::AbstractMatrix,
                        η::AbstractMatrix,
                        G::Grid,
                        Diag::DiagnosticVars,
                        Forc::Forcing)

    @unpack h,h_u,h_v,U,V,dUdx,dVdy = Diag.VolumeFluxes
    @unpack g = C
    @unpack ep = G

    # mass or volume flux U,V = uH,vH; h_u, h_v are actually H_u, H_v
    Uflux!(U,u,h_u,ep)
    Vflux!(V,v,h_v)

    # divergence of mass flux
    ∂x!(dUdx,U)
    ∂y!(dVdy,V)

    # Pressure gradient
    @unpack dpdx,dpdy = Diag.Bernoulli
    ∂x!(dpdx,g*η)
    ∂y!(dpdy,g*η)

    # Coriolis force
    @unpack qhv,qhu,v_u,u_v = Diag.Vorticity
    @unpack f_u,f_v = G
    Ixy!(v_u,v)
    Ixy!(u_v,u)
    fv!(qhv,f_u,v_u)
    fu!(qhu,f_v,u_v)

    # adding the terms
    @unpack du,dv,dη = Diag.Tendencies
    @unpack Fx,Fy = Forc
    momentum_u!(du,qhv,dpdx,Fx)
    momentum_v!(dv,qhu,dpdy,Fy)
    continuity!(dη,dUdx,dVdy)
end

""" Update advective and Coriolis tendencies."""
function advection_coriolis!(   u::AbstractMatrix,
                                v::AbstractMatrix,
                                η::AbstractMatrix,
                                P::Parameter,
                                G::Grid,
                                Diag::DiagnosticVars,
                                Forc::Forcing)

    @unpack h = Diag.VolumeFluxes
    @unpack H = Forc
    @unpack q,h_q,dvdx,dudy = Diag.Vorticity
    @unpack u²,v²,KEu,KEv = Diag.Bernoulli
    @unpack ep,f_q = G

    if G.nstep_advcor > 0
        thickness!(h,η,H)
    end

    Ixy!(h_q,h)

    # off-diagonals of stress tensor ∇(u,v)
    ∂x!(dvdx,v)
    ∂y!(dudy,u)

    # non-linear part of the Bernoulli potential
    speed!(u²,v²,u,v)
    Ix!(KEu,u²)
    Iy!(KEv,v²)

    # Potential vorticity update
    PV!(G,Diag)
    
    @unpack q = Diag.Vorticity
    # Linear combinations of the potential vorticity q
    if P.adv_scheme == "Sadourny"
        @unpack q_u,q_v = Diag.Vorticity
        Iy!(q_u,q)
        Ix!(q_v,q)
    elseif P.adv_scheme == "ArakawaHsu"
        @unpack qα,qβ,qγ,qδ = Diag.ArakawaHsu
        AHα!(qα,q)
        AHβ!(qβ,q)
        AHγ!(qγ,q)
        AHδ!(qδ,q)
    end
end

"""Layer thickness h obtained by adding sea surface height η to bottom height H."""
function thickness!(h::AbstractMatrix,η::AbstractMatrix,H::AbstractMatrix)
    m,n = size(h)
    @boundscheck (m,n) == size(η) || throw(BoundsError())
    @boundscheck (m,n) == size(H) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            h[i,j] = η[i,j] + H[i,j]
        end
    end
end

"""Zonal mass flux U = uh."""
function Uflux!(U::AbstractMatrix,
                u::AbstractMatrix,
                h_u::AbstractMatrix,
                ep::Int)

    m,n = size(U)
    @boundscheck (m,n) == size(h_u) || throw(BoundsError())
    @boundscheck (m+2+ep,n+2) == size(u) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            U[i,j] = u[1+ep+i,1+j]*h_u[i,j]
        end
    end
end

"""Meridional mass flux V = vh."""
function Vflux!(V::AbstractMatrix,v::AbstractMatrix,h_v::AbstractMatrix)
    m,n = size(V)
    @boundscheck (m,n) == size(h_v) || throw(BoundsError())
    @boundscheck (m+2,n+2) == size(v) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            V[i,j] = v[i+1,j+1]*h_v[i,j]
        end
    end
end

"""Squared velocities u²,v²."""
function speed!(u²::AbstractMatrix,
                v²::AbstractMatrix,
                u::AbstractMatrix,
                v::AbstractMatrix)
    m,n = size(u²)
    @boundscheck (m,n) == size(u) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            u²[i,j] = u[i,j]^2
        end
    end

    m,n = size(v²)
    @boundscheck (m,n) == size(v) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            v²[i,j] = v[i,j]^2
        end
    end
end

"""Bernoulli potential p = 1/2*(u² + v²) + gη."""
function bernoulli!(p::Array{T,2},
                    KEu::Array{T,2},
                    KEv::Array{T,2},
                    η::Array{T,2},
                    g::T,
                    ep::Int) where {T<:AbstractFloat}
    m,n = size(p)
    @boundscheck (m+ep,n+2) == size(KEu) || throw(BoundsError())
    @boundscheck (m+2,n) == size(KEv) || throw(BoundsError())
    @boundscheck (m,n) == size(η) || throw(BoundsError())

    one_half = T(0.5)

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            p[i,j] = one_half*(KEu[i+ep,j+1] + KEv[i+1,j]) + g*η[i,j]
        end
    end
end

"""Coriolis term f*v. """
function fv!(   qhv::AbstractMatrix,
                f_u::AbstractMatrix,
                v_u::AbstractMatrix,
                ep::Int)

    m,n = size(qhv)
    @boundscheck (m,n) == size(f_u) || throw(BoundsError())
    @boundscheck (m+4-ep,n+2) == size(v_u) || throw(BoundsError())
    @boundscheck ep == 1 || ep == 0 || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            qhv[i,j] = f_u[i,j]*v_u[i+2-ep,j+1]
        end
    end
end

"""Coriolis term f*u. """
function fu!(   qhu::AbstractMatrix,
                f_v::AbstractMatrix,
                u_v::AbstractMatrix,
                ep::Int)

    m,n = size(qhu)
    @boundscheck (m,n) == size(f_v) || throw(BoundsError())
    @boundscheck (m+2+ep,n+4) == size(u_v) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            qhu[i,j] = f_v[i,j]*u_v[i+1+ep,j+2]
        end
    end
end

"""Sum up the tendencies of the non-diffusive right-hand side for the u-component."""
function momentum_u!(   du::AbstractMatrix,
                        qhv::AbstractMatrix,
                        dpdx::AbstractMatrix,
                        Fx::AbstractMatrix,
                        ep::Int)

    m,n = size(du) .- (4,4)     # cut off the halo
    @boundscheck (m,n) == size(qhv) || throw(BoundsError())
    @boundscheck (m+2-ep,n+2) == size(dpdx) || throw(BoundsError())
    @boundscheck (m,n) == size(Fx) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            du[i+2,j+2] = qhv[i,j] - dpdx[i+1-ep,j+1] + Fx[i,j]
        end
    end
end

"""Sum up the tendencies of the non-diffusive right-hand side for the v-component."""
function momentum_v!(   dv::AbstractMatrix,
                        qhu::AbstractMatrix,
                        dpdy::AbstractMatrix,
                        Fy::AbstractMatrix)

    m,n = size(dv) .- (4,4)     # cut off the halo
    @boundscheck (m,n) == size(qhu) || throw(BoundsError())
    @boundscheck (m+2,n+2) == size(dpdy) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
             dv[i+2,j+2] = -qhu[i,j] - dpdy[i+1,j+1] + Fy[i,j]
        end
    end
end

"""Tendencies du,dv,dη of the non-diffusive right-hand side."""
function rhs!(du,dv,dη,u,v,η,Fx,f_q,H,η_ref,
            dvdx,dudy,dpdx,dpdy,
            p,u²,v²,KEu,KEv,dUdx,dVdy,
            h,h_u,h_v,h_q,U,V,U_v,V_u,
            qhv,qhu,q,q_u,q_v,
            qα,qβ,qγ,qδ)

    # layer thickness
    thickness!(h,η,H)
    Ix!(h_u,h)
    Iy!(h_v,h)
    Ixy!(h_q,h)

    # mass or volume flux U,V = uh,vh
    Uflux!(U,u,h_u)
    Vflux!(V,v,h_v)

    # off-diagonals of stress tensor ∇(u,v), ∇(U,V)
    ∂x!(dvdx,v)
    ∂y!(dudy,u)
    ∂x!(dUdx,U)
    ∂y!(dVdy,V)

    # Bernoulli potential
    speed!(u²,v²,u,v)
    Ix!(KEu,u²)
    Iy!(KEv,v²)
    Bernoulli!(p,KEu,KEv,η)
    ∂x!(dpdx,p)
    ∂y!(dpdy,p)

    # Potential vorticity and advection thereof
    PV!(q,f_q,dvdx,dudy,h_q)
    PV_adv!(qhv,qhu,q,qα,qβ,qγ,qδ,q_u,q_v,U,V,V_u,U_v)

    # adding the terms
    momentum_u!(du,qhv,dpdx,Fx)
    momentum_v!(dv,qhu,dpdy)
    continuity!(dη,dUdx,dVdy,η,η_ref)
end

"""Layer thickness h obtained by adding sea surface height η to bottom height H."""
function thickness!(h,η,H)
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
function Uflux!(U,u,h_u)
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
function Vflux!(V,v,h_v)
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
function speed!(u²,v²,u,v)
    # u squared
    m,n = size(u²)
    @boundscheck (m,n) == size(u) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            u²[i,j] = u[i,j]^2
        end
    end

    # v squared
    m,n = size(v²)
    @boundscheck (m,n) == size(v) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            v²[i,j] = v[i,j]^2
        end
    end
end

"""Bernoulli potential p = 1/2*(u² + v²) + gη."""
function Bernoulli!(p,KEu,KEv,η)
    m,n = size(p)
    @boundscheck (m+ep,n+2) == size(KEu) || throw(BoundsError())
    @boundscheck (m+2,n) == size(KEv) || throw(BoundsError())
    @boundscheck (m,n) == size(η) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            p[i,j] = one_half*(KEu[i+ep,j+1] + KEv[i+1,j]) + g*η[i,j]
        end
    end
end

"""Sum up the tendencies of the non-diffusive right-hand side for the u-component."""
function momentum_u!(du,qhv,dpdx,Fx)
    m,n = size(du) .- (2*halo,2*halo) # cut off the halo
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
function momentum_v!(dv,qhu,dpdy)
    m,n = size(dv) .- (2*halo,2*halo) # cut off the halo
    @boundscheck (m,n) == size(qhu) || throw(BoundsError())
    @boundscheck (m+2,n+2) == size(dpdy) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
             dv[i+2,j+2] = -qhu[i,j] - dpdy[i+1,j+1]
        end
    end
end

"""Continuity equation's right-hand side with surface relaxation
-∂x(uh) - ∂y(vh) + γ*(η_ref - η)."""
function continuity_surf_forc!(dη,dUdx,dVdy,η,η_ref)
    m,n = size(dη) .- (2*haloη,2*haloη)
    @boundscheck (m,n+2) == size(dUdx) || throw(BoundsError())
    @boundscheck (m+2,n) == size(dVdy) || throw(BoundsError())
    @boundscheck (m+2,n+2) == size(η) || throw(BoundsError())
    @boundscheck (m,n) == size(η_ref) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            dη[i+1,j+1] = -(dUdx[i,j+1] + dVdy[i+1,j]) + γ*(η_ref[i,j]-η[i+1,j+1])
        end
    end
end

"""Continuity equation's right-hand side -∂x(uh) - ∂y(vh) without forcing."""
function continuity_itself!(dη,dUdx,dVdy,η,η_ref)
    m,n = size(dη) .- (2*haloη,2*haloη)
    @boundscheck (m,n+2) == size(dUdx) || throw(BoundsError())
    @boundscheck (m+2,n) == size(dVdy) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            dη[i+1,j+1] = -(dUdx[i,j+1] + dVdy[i+1,j])
        end
    end
end

if surface_forcing
    continuity! = continuity_surf_forc!
else
    continuity! = continuity_itself!
end

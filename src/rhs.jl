function rhs!(du,dv,dη,u,v,η,Fx,f_q,H,
            dudx,dvdy,dvdx,dudy,dpdx,dpdy,
            p,u²,v²,KEu,KEv,dUdx,dVdy,
            h,h_u,h_v,h_q,U,V,U_v,V_u,
            qhv,qhu,q,q_u,q_v,
            qα,qβ,qγ,qδ,
            sqrtKE,sqrtKE_u,sqrtKE_v,Bu,Bv,
            DS,DS_q,DT,νSmag,νSmag_q,
            Lu,Lv,dLudx,dLudy,dLvdx,dLvdy,
            S11,S12,S21,S22,
            LLu1,LLu2,LLv1,LLv2)

    # compute the tendencies du,dv,dη of the right-hand side

    # layer thickness
    thickness!(h,η,H)
    Ix!(h_u,h)
    Iy!(h_v,h)
    Ixy!(h_q,h)

    # mass or volume flux U,V = uh,vh
    Uflux!(U,u,h_u)
    Vflux!(V,v,h_v)

    # stress tensor ∇u⃗
    ∂x!(dudx,u)
    ∂y!(dvdy,v)
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

    # Potential vorticity
    PV!(q,f_q,dvdx,dudy,h_q)

    # Sadourny, 1975 enstrophy conserving scheme
    if adv_scheme == "Sadourny"
        PV_adv_Sadourny!(qhv,qhu,q,q_u,q_v,U,V,V_u,U_v)
    elseif adv_scheme == "ArakawaHsu"
        PV_adv_ArakawaHsu!(qhv,qhu,q,qα,qβ,qγ,qδ,U,V)
    end

    # bottom drag
    bottom_drag!(Bu,Bv,KEu,KEv,sqrtKE,sqrtKE_u,sqrtKE_v,u,v,h_u,h_v)

    # Smagorinsky-like biharmonic diffusion
    Smagorinsky_coeff!(νSmag,νSmag_q,DS,DS_q,DT,dudx,dvdy,dudy,dvdx)
    stress_tensor!(dLudx,dLudy,dLvdx,dLvdy,Lu,Lv,u,v)
    viscous_tensor!(S11,S12,S21,S22,νSmag,νSmag_q,dLudx,dLudy,dLvdx,dLvdy)

    ∂x!(LLu1,S11)
    ∂y!(LLu2,S12)
    ∂x!(LLv1,S21)
    ∂y!(LLv2,S22)

    # adding the terms
    momentum_u!(du,qhv,dpdx,Bu,LLu1,LLu2,Fx)
    momentum_v!(dv,qhu,dpdy,Bv,LLv1,LLv2)
    continuity!(dη,dUdx,dVdy)
end

function thickness!(h,η,H)
    # h = η + H
    m,n = size(h)
    @boundscheck (m,n) == size(η) || throw(BoundsError())
    @boundscheck (m,n) == size(H) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            h[i,j] = η[i,j] + H[i,j]
        end
    end
end

function Uflux!(U,u,h_u)
    # U = uh
    m,n = size(U)
    @boundscheck (m,n) == size(h_u) || throw(BoundsError())
    @boundscheck (m+2+ep,n+2) == size(u) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            U[i,j] = u[1+ep+i,1+j]*h_u[i,j]
        end
    end
end

function Vflux!(V,v,h_v)
    # V = vh
    m,n = size(V)
    @boundscheck (m,n) == size(h_v) || throw(BoundsError())
    @boundscheck (m+2,n+2) == size(v) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            V[i,j] = v[i+1,j+1]*h_v[i,j]
        end
    end
end

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

function Bernoulli!(p,KEu,KEv,η)
    # Bernoulli potential p = 1/2*(u² + v²) + gη
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

function PV!(q,f_q,dvdx,dudy,h_q)
    # Potential vorticity q = (f + ∂v/∂x - ∂u/∂y)/h
    m,n = size(q)
    @boundscheck (m,n) == size(f_q) || throw(BoundsError())
    @boundscheck (m+2,n+2) == size(dvdx) || throw(BoundsError())
    @boundscheck (m+2+ep,n+2) == size(dudy) || throw(BoundsError())
    @boundscheck (m,n) == size(h_q) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            q[i,j] = (f_q[i,j] + dvdx[i+1,j+1] - dudy[i+1+ep,j+1]) / h_q[i,j]
        end
    end
end

function PV_adv_Sadourny!(qhv,qhu,q,q_u,q_v,U,V,V_u,U_v)
    #= Advection of potential voriticity qhv,qhu as in Sadourny, 1975
    enstrophy conserving scheme =#

    Iy!(q_u,q)
    Ix!(q_v,q)
    Ixy!(V_u,V)
    Ixy!(U_v,U)

    @views qhv .= q_u[2-ep:end-1,:].*V_u[2-ep:end-1,:]
    @views qhu .= q_v[:,2:end-1].*U_v[:,2:end-1]
end

function PV_adv_ArakawaHsu!(qhv,qhu,q,qα,qβ,qγ,qδ,U,V)
    #= Advection of potential vorticity qhv,qhu as in Arakawa and Hsu, 1990
    Energy and enstrophy conserving (in the limit of non-divergent mass flux) scheme with τ = 0. =#

    # Linear combinations α,β,γ,δ of potential vorticity q
    AHα!(qα,q)
    AHβ!(qβ,q)
    AHγ!(qγ,q)
    AHδ!(qδ,q)

    # Linear combinations of q and V=hv to yield qhv
    m,n = size(qhv)
    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            qhv[i,j] = qα[1-ep+i,j].*V[2-ep+i,j+1] .+ qβ[1-ep+i,j].*V[1-ep+i,j+1] .+ qγ[1-ep+i,j].*V[1-ep+i,j] .+ qδ[1-ep+i,j].*V[2-ep+i,j]
        end
    end

    # Linear combinations of q and U=hu to yield qhu
    m,n = size(qhu)
    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            qhu[i,j] = qα[i,j].*U[i,j+1] .+ qβ[i+1,j].*U[i+1,j+1] .+ qγ[i+1,j+1].*U[i+1,j+2] .+ qδ[i,j+1].*U[i,j+2]
        end
    end
end

function bottom_drag!(Bu,Bv,KEu,KEv,sqrtKE,sqrtKE_u,sqrtKE_v,u,v,h_u,h_v)
    # quadratic bottom drag Bu,Bv = c_D/h * | u⃗ | * u⃗

    # sqrt of KE, which is actually the kinetic energy
    m,n = size(sqrtKE)
    @boundscheck (m+ep,n+2) == size(KEu) || throw(BoundsError())
    @boundscheck (m+2,n) == size(KEv) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            sqrtKE[i,j] = sqrt(KEu[i+ep,j+1] + KEv[i+1,j])
        end
    end

    Ix!(sqrtKE_u,sqrtKE)
    Iy!(sqrtKE_v,sqrtKE)

    m,n = size(Bu)
    @boundscheck (m,n) == size(sqrtKE_u) || throw(BoundsError())
    @boundscheck (m,n) == size(h_u) || throw(BoundsError())
    @boundscheck (m+2+ep,n+2) == size(u) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            Bu[i,j] = c_D*sqrtKE_u[i,j] * u[i+1+ep,j+1] / h_u[i,j]
        end
    end

    m,n = size(Bv)
    @boundscheck (m,n) == size(sqrtKE_v) || throw(BoundsError())
    @boundscheck (m,n) == size(h_v) || throw(BoundsError())
    @boundscheck (m+2,n+2) == size(v) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            Bv[i,j] = c_D*sqrtKE_v[i,j] * v[i+1,j+1] / h_v[i,j]
        end
    end
end


function Smagorinsky_coeff!(νSmag,νSmag_q,DS,DS_q,DT,dudx,dvdy,dudy,dvdx)
    # νSmag = cSmag * |D|, where deformation rate |D| = √((∂u/∂x - ∂v/∂y)^2 + (∂u/∂y + ∂v/∂x)^2)
    # the grid spacing Δ is omitted here as the operators are dimensionless

    # horizontal shearing strain squared
    @views DS_q .= (dudy[1+ep:end,:] .+ dvdx).^2
    Ixy!(DS,DS_q)

    # horizontal tension squared
    @views DT .= (dudx[1+ep:end,2:end-1] + dvdy[2:end-1,:]).^2

    # Smagorinsky coefficient times deformation rate
    @views νSmag .= cSmag*sqrt.(DS .+ DT)
    Ixy!(νSmag_q,νSmag)
end

function stress_tensor!(dLudx,dLudy,dLvdx,dLvdy,Lu,Lv,u,v)
    # Biharmonic stress tensor ∇∇²u⃗ = (∂/∂x(∇²u), ∂/∂y(∇²u); ∂/∂x(∇²v), ∂/∂y(∇²v))
    ∇²!(Lu,u)
    ∇²!(Lv,v)
    ∂x!(dLudx,Lu)
    ∂y!(dLudy,Lu)
    ∂x!(dLvdx,Lv)
    ∂y!(dLvdy,Lv)
end

function viscous_tensor!(S11,S12,S21,S22,νSmag,νSmag_q,dLudx,dLudy,dLvdx,dLvdy)
    # Biharmonic stress tensor times Smagorinsky coefficient
    # νSmag * ∇∇² ⃗u = (S11, S12; S21, S22)
    m,n = size(S11)
    @boundscheck (m+2-ep,n) == size(νSmag) || throw(BoundsError())
    @boundscheck (m,n) == size(dLudx) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            S11[i,j] = νSmag[i+1-ep,j] * dLudx[i,j]
        end
    end

    m,n = size(S12)
    @boundscheck (m,n) == size(νSmag_q) || throw(BoundsError())
    @boundscheck (m+ep,n) == size(dLudy) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            S12[i,j] = νSmag_q[i,j] * dLudy[i+ep,j]
        end
    end

    m,n = size(S21)
    @boundscheck (m,n) == size(νSmag_q) || throw(BoundsError())
    @boundscheck (m,n) == size(dLvdx) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            S21[i,j] = νSmag_q[i,j] * dLvdx[i,j]
        end
    end

    m,n = size(S22)
    @boundscheck (m,n+2) == size(νSmag) || throw(BoundsError())
    @boundscheck (m,n) == size(dLvdy) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            S22[i,j] = νSmag[i,j+1] * dLvdy[i,j]
        end
    end
end

function momentum_u!(du,qhv,dpdx,Bu,LLu1,LLu2,Fx)
    # Sum up the tendency terms of the right-hand side for the u-component
    m,n = size(du) .- (2*halo,2*halo) # cut off the halo
    @boundscheck (m,n) == size(qhv) || throw(BoundsError())
    @boundscheck (m+2-ep,n+2) == size(dpdx) || throw(BoundsError())
    @boundscheck (m+2-ep,n+2) == size(Bu) || throw(BoundsError())
    @boundscheck (m,n+2) == size(LLu1) || throw(BoundsError())
    @boundscheck (m+2-ep,n) == size(LLu2) || throw(BoundsError())
    @boundscheck (m,n) == size(Fx) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            du[i+2,j+2] = qhv[i,j] - dpdx[i+1-ep,j+1] - Bu[i+1-ep,j+1] + LLu1[i,j+1] + LLu2[i+1-ep,j] + Fx[i,j]
        end
    end
end

function momentum_v!(dv,qhu,dpdy,Bv,LLv1,LLv2)
    # Sum up the tendency terms of the right-hand side for the v-component
    m,n = size(dv) .- (2*halo,2*halo) # cut off the halo
    @boundscheck (m,n) == size(qhu) || throw(BoundsError())
    @boundscheck (m+2,n+2) == size(dpdy) || throw(BoundsError())
    @boundscheck (m+2,n+2) == size(Bv) || throw(BoundsError())
    @boundscheck (m,n+2) == size(LLv1) || throw(BoundsError())
    @boundscheck (m+2,n) == size(LLv2) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
             dv[i+2,j+2] = -qhu[i,j] - dpdy[i+1,j+1] - Bv[i+1,j+1] + LLv1[i,j+1] + LLv2[i+1,j]
        end
    end
end

function continuity!(dη,dUdx,dVdy)
    # Continuity equation's right-hand side -∂x(uh) - ∂y(vh)
    m,n = size(dη) .- (2*haloη,2*haloη)
    @boundscheck (m,n+2) == size(dUdx) || throw(BoundsError())
    @boundscheck (m+2,n) == size(dVdy) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            dη[i+1,j+1] = -(dUdx[i,j+1] + dVdy[i+1,j])
        end
    end
end

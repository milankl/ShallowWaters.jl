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

    @views h .= η .+ H
    Ix!(h_u,h)
    Iy!(h_v,h)
    Ixy!(h_q,h)

    Uflux!(U,u,h_u)
    Vflux!(V,v,h_v)

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

function Uflux!(U,u,h_u)
    # U = uh
    @views U .= u[2+ep:end-1,2:end-1].*h_u
end

function Vflux!(V,v,h_v)
    # V = vh
    @views V .= v[2:end-1,2:end-1].*h_v
end

function speed!(u²,v²,u,v)
    @views u² .= u.^2
    @views v² .= v.^2
end

function Bernoulli!(p,KEu,KEv,η)
    # Bernoulli potential p = 1/2*(u² + v²) + gη
    @views p .= one_half*(KEu[1+ep:end,2:end-1] .+ KEv[2:end-1,:]) .+ g*η
end

function PV!(q,f_q,dvdx,dudy,h_q)
    # Potential vorticity q = (f + ∂v/∂x - ∂u/∂y)/h
    @views q .= (f_q .+ dvdx[2:end-1,2:end-1] .- dudy[2+ep:end-1,2:end-1]) ./ h_q
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

    @views qhv .= qα[2-ep:end,:].*V[3-ep:end-1,2:end]
    @views qhv .+= qβ[2-ep:end-1,:].*V[2-ep:end-2,2:end]
    @views qhv .+= qγ[2-ep:end-1,:].*V[2-ep:end-2,1:end-1]
    @views qhv .+= qδ[2-ep:end,:].*V[3-ep:end-1,1:end-1]

    @views qhu .= qα[:,1:end-1].*U[1:end-1,2:end-2]
    @views qhu .+= qβ[2:end,1:end-1].*U[2:end,2:end-2]
    @views qhu .+= qγ[2:end,2:end].*U[2:end,3:end-1]
    @views qhu .+= qδ[:,2:end].*U[1:end-1,3:end-1]

end

function bottom_drag!(Bu,Bv,KEu,KEv,sqrtKE,sqrtKE_u,sqrtKE_v,u,v,h_u,h_v)
    # quadratic bottom drag Bu,Bv = c_D/h * | u⃗ | * u⃗
    @views sqrtKE .= sqrt.(KEu[1+ep:end,2:end-1] .+ KEv[2:end-1,:])
    Ix!(sqrtKE_u,sqrtKE)
    Iy!(sqrtKE_v,sqrtKE)
    @views Bu .= c_D*sqrtKE_u .* u[2+ep:end-1,2:end-1] ./ h_u
    @views Bv .= c_D*sqrtKE_v .* v[2:end-1,2:end-1] ./ h_v
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
    @views S11 .= νSmag[2-ep:end-1,:].*dLudx
    @views S12 .= νSmag_q.*dLudy[1+ep:end,:]
    @views S21 .= νSmag_q.*dLvdx
    @views S22 .= νSmag[:,2:end-1].*dLvdy
end

function momentum_u!(du,qhv,dpdx,Bu,LLu1,LLu2,Fx)
    # Sum up the tendency terms of the right-hand side for the u-component
    @views du[3:end-2,3:end-2] .= qhv .- dpdx[2-ep:end-1,2:end-1] .- Bu[2-ep:end-1,2:end-1] .+ LLu1[:,2:end-1] .+ LLu2[2-ep:end-1,:] .+ Fx
end

function momentum_v!(dv,qhu,dpdy,Bv,LLv1,LLv2)
    # Sum up the tendency terms of the right-hand side for the v-component
    @views dv[3:end-2,3:end-2] .= -qhu .- dpdy[2:end-1,2:end-1] .- Bv[2:end-1,2:end-1] .+ LLv1[:,2:end-1] + LLv2[2:end-1,:]
end

function continuity!(dη,dUdx,dVdy)
    # Continuity equation's right-hand side -∂x(uh) - ∂y(vh)
    m,n = size(dη)

    @inbounds for i ∈ 2:n-1
        for j ∈ 2:m-1
            dη[j,i] = -(dUdx[j-1,i] + dVdy[j,i-1])
        end
    end
end

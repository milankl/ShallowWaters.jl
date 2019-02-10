"""Biharmonic diffusion operator with constant viscosity coefficient
    Viscosity = ν∇⁴ ⃗u. Although constant, the coefficient is actually inside
    Viscosity = ∇⋅ν∇∇² ⃗u."""
function diffusion_constant!(dudx,dudy,dvdx,dvdy,DS,DS_q,DT,νSmag,νSmag_q,Lu,Lv,dLudx,dLudy,dLvdx,dLvdy,
        S11,S12,S21,S22,LLu1,LLu2,LLv1,LLv2,u,v)
    stress_tensor!(dLudx,dLudy,dLvdx,dLvdy,Lu,Lv,u,v)
    viscous_tensor_Constant!(S11,S12,S21,S22,dLudx,dLudy,dLvdx,dLvdy)

    ∂x!(LLu1,S11)
    ∂y!(LLu2,S12)
    ∂x!(LLv1,S21)
    ∂y!(LLv2,S22)
end

""" Smagorinsky-like biharmonic viscosity
    Viscosity = ∇ ⋅ (cSmag Δ⁴ |D| ∇∇² ⃗u)
The Δ⁴-scaling is omitted as gradient operators are dimensionless."""
function diffusion_Smagorinsky!(dudx,dudy,dvdx,dvdy,DS,DS_q,DT,νSmag,νSmag_q,Lu,Lv,dLudx,dLudy,dLvdx,dLvdy,
        S11,S12,S21,S22,LLu1,LLu2,LLv1,LLv2,u,v)
    ∂x!(dudx,u)
    ∂y!(dvdy,v)
    ∂x!(dvdx,v)
    ∂y!(dudy,u)

    # biharmonic diffusion
    stress_tensor!(dLudx,dLudy,dLvdx,dLvdy,Lu,Lv,u,v)
    Smagorinsky_coeff!(νSmag,νSmag_q,DS,DS_q,DT,dudx,dvdy,dudy,dvdx)
    viscous_tensor_Smagorinsky!(S11,S12,S21,S22,νSmag,νSmag_q,dLudx,dLudy,dLvdx,dLvdy)

    ∂x!(LLu1,S11)
    ∂y!(LLu2,S12)
    ∂x!(LLv1,S21)
    ∂y!(LLv2,S22)
end

"""νSmag = cSmag * |D|, where deformation rate |D| = √((∂u/∂x - ∂v/∂y)^2 + (∂u/∂y + ∂v/∂x)^2).
The grid spacing Δ is omitted here as the operators are dimensionless."""
function Smagorinsky_coeff!(νSmag,νSmag_q,DS,DS_q,DT,dudx,dvdy,dudy,dvdx)
    # horizontal shearing strain squared
    m,n = size(DS_q)
    @boundscheck (m+ep,n) == size(dudy) || throw(BoundsError())
    @boundscheck (m,n) == size(dvdx) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            DS_q[i,j] = (dudy[i+ep,j] + dvdx[i,j])^2
        end
    end

    Ixy!(DS,DS_q)

    # horizontal tension squared
    m,n = size(DT)
    @boundscheck (m+ep,n+2) == size(dudx) || throw(BoundsError())
    @boundscheck (m+2,n) == size(dvdy) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            DT[i,j] = (dudx[i+ep,j+1] + dvdy[i+1,j])^2
        end
    end

    # viscosity = Smagorinsky coefficient times deformation rate
    m,n = size(νSmag)
    @boundscheck (m,n) == size(DS) || throw(BoundsError())
    @boundscheck (m,n) == size(DT) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            νSmag[i,j] = cSmag*sqrt(DS[i,j] + DT[i,j])
        end
    end

    Ixy!(νSmag_q,νSmag)
end

"""Biharmonic stress tensor ∇∇²(u,v) = (∂/∂x(∇²u), ∂/∂y(∇²u); ∂/∂x(∇²v), ∂/∂y(∇²v))"""
function stress_tensor!(dLudx,dLudy,dLvdx,dLvdy,Lu,Lv,u,v)
    ∇²!(Lu,u)
    ∇²!(Lv,v)
    ∂x!(dLudx,Lu)
    ∂y!(dLudy,Lu)
    ∂x!(dLvdx,Lv)
    ∂y!(dLvdy,Lv)
end

"""Biharmonic stress tensor times Smagorinsky coefficient
νSmag * ∇∇² ⃗u = (S11, S12; S21, S22)."""
function viscous_tensor_Smagorinsky!(S11,S12,S21,S22,νSmag,νSmag_q,dLudx,dLudy,dLvdx,dLvdy)
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

"""Biharmonic stress tensor times constant viscosity coefficient
νB * ∇∇² ⃗u = (S11, S12; S21, S22)"""
function viscous_tensor_Constant!(S11,S12,S21,S22,dLudx,dLudy,dLvdx,dLvdy)
    m,n = size(S11)
    @boundscheck (m,n) == size(dLudx) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            S11[i,j] = νB * dLudx[i,j]
        end
    end

    m,n = size(S12)
    @boundscheck (m+ep,n) == size(dLudy) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            S12[i,j] = νB * dLudy[i+ep,j]
        end
    end

    m,n = size(S21)
    @boundscheck (m,n) == size(dLvdx) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            S21[i,j] = νB * dLvdx[i,j]
        end
    end

    m,n = size(S22)
    @boundscheck (m,n) == size(dLvdy) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            S22[i,j] = νB * dLvdy[i,j]
        end
    end
end

"""Update u with bottom friction tendency (Bu,Bv) and biharmonic viscosity."""
function add_drag_diff_tendencies!(u,v,Bu,Bv,LLu1,LLu2,LLv1,LLv2)
    m,n = size(u) .- (2*halo,2*halo)
    @boundscheck (m+2-ep,n+2) == size(Bu) || throw(BoundsError())
    @boundscheck (m,n+2) == size(LLu1) || throw(BoundsError())
    @boundscheck (m+2-ep,n) == size(LLu2) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            u[i+2,j+2] += nΔt_diff*(Bu[i+1-ep,j+1] + LLu1[i,j+1] + LLu2[i+1-ep,j])
        end
    end

    m,n = size(v) .- (2*halo,2*halo) # cut off the halo
    @boundscheck (m+2,n+2) == size(Bv) || throw(BoundsError())
    @boundscheck (m,n+2) == size(LLv1) || throw(BoundsError())
    @boundscheck (m+2,n) == size(LLv2) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
             v[i+2,j+2] += nΔt_diff*(Bv[i+1,j+1] + LLv1[i,j+1] + LLv2[i+1,j])
        end
    end
end

if diffusion == "Constant"
    diffusive! = diffusion_constant!
elseif diffusion == "Smagorinsky"
    diffusive! = diffusion_Smagorinsky!
else
    throw(error("Diffusion not correctly specified. Only Constant or Smagorinsky allowed."))
end

function bottom_drag_quad!(Bu,Bv,KEu,KEv,sqrtKE,sqrtKE_u,sqrtKE_v,u,v,η,H,h,u²,v²,h_u,h_v)
    # quadratic bottom drag Bu,Bv = c_D/h * | u⃗ | * u⃗

    thickness!(h,η,H)
    Ix!(h_u,h)
    Iy!(h_v,h)

    speed!(u²,v²,u,v)
    Ix!(KEu,u²)
    Iy!(KEv,v²)

    # sqrt of KE, which is actually the kinetic energy without the 0.5 factor
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

"""Linear bottom drag computed as r_B*(u,v). r_B is negative and contains the
grid spacing Δ as gradient operators are dimensionless."""
function bottom_drag_lin!(Bu,Bv,KEu,KEv,sqrtKE,sqrtKE_u,sqrtKE_v,u,v,η,H,h,u²,v²,h_u,h_v)
    m,n = size(Bu)
    @boundscheck (m+2+ep,n+2) == size(u) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            Bu[i,j] = r_B * u[i+1+ep,j+1]
        end
    end

    m,n = size(Bv)
    @boundscheck (m+2,n+2) == size(v) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            Bv[i,j] = r_B * v[i+1,j+1]
        end
    end
end

function no_bottom_drag!(Bu,Bv,KEu,KEv,sqrtKE,sqrtKE_u,sqrtKE_v,u,v,η,H,h,u²,v²,h_u,h_v)
    Bu .= zeero
    Bv .= zeero
end

if bottom_friction == "linear"
    bottom_drag! = bottom_drag_lin!
elseif bottom_friction == "quadratic"
    bottom_drag! = bottom_drag_quad!
elseif bottom_friction == "none"
    bottom_drag! = no_bottom_drag!
else
    throw(error("Bottom friction not correctly specified. Only linear, quadratic or none allowed."))
end

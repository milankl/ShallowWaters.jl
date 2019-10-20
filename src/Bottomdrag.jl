"""Transit function to call the specified bottom drag function."""
function bottom_drag!(  u::AbstractMatrix,
                        v::AbstractMatrix,
                        η::AbstractMatrix,
                        Diag::DiagnosticVars,
                        S::ModelSetup)

    @unpack bottom_drag = S.parameters

    if bottom_drag == "quadratic"
        bottom_drag_quadratic!(u,v,η,Diag,S)
    elseif bottom_drag == "linear"
        bottom_drag_linear!(u,v,η,Diag,S)
    elseif bottom_drag == "none"
        no_bottom_drag!(u,v,η,Diag,S)
    end
end

"""Quadratic bottom drag Bu,Bv = cD/h * | u⃗ | * u⃗"""
function bottom_drag_quadratic!(u::AbstractMatrix,
                                v::AbstractMatrix,
                                η::AbstractMatrix,
                                Diag::DiagnosticVars,
                                S::ModelSetup)

    @unpack h,h_u,h_v = Diag.VolumeFluxes
    @unpack u²,v²,KEu,KEv = Diag.Bernoulli
    @unpack Bu,Bv,sqrtKE,sqrtKE_u,sqrtKE_v = Diag.Bottomdrag
    @unpack ep = S.grid
    @unpack cD = S.constants
    @unpack H = S.forcing

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
            Bu[i,j] = cD*sqrtKE_u[i,j] * u[i+1+ep,j+1] / h_u[i,j]
        end
    end

    m,n = size(Bv)
    @boundscheck (m,n) == size(sqrtKE_v) || throw(BoundsError())
    @boundscheck (m,n) == size(h_v) || throw(BoundsError())
    @boundscheck (m+2,n+2) == size(v) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            Bv[i,j] = cD*sqrtKE_v[i,j] * v[i+1,j+1] / h_v[i,j]
        end
    end
end

"""Linear bottom drag computed as r_B*(u,v). r_B is negative and contains the
grid spacing Δ as gradient operators are dimensionless."""
function bottom_drag_linear!(   u::AbstractMatrix,
                                v::AbstractMatrix,
                                η::AbstractMatrix,
                                Diag::DiagnosticVars,
                                S::ModelSetup)

    @unpack Bu,Bv = Diag.Bottomdrag
    @unpack ep = S.grid
    @unpack rD = S.constants

    m,n = size(Bu)
    @boundscheck (m+2+ep,n+2) == size(u) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            Bu[i,j] = rD * u[i+1+ep,j+1]
        end
    end

    m,n = size(Bv)
    @boundscheck (m+2,n+2) == size(v) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            Bv[i,j] = rD * v[i+1,j+1]
        end
    end
end

"""No bottom drag."""
function no_bottom_drag!(   u::AbstractMatrix,
                            v::AbstractMatrix,
                            η::AbstractMatrix,
                            Diag::DiagnosticVars,
                            S::ModelSetup)

    @unpack Bu,Bv = Diag.Bottomdrag
    Bu .= zero(eltype(Bu))
    Bv .= zero(eltype(Bu))
end

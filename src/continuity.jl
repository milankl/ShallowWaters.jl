"""Continuity equation's right-hand side with surface relaxation
-∂x(uh) - ∂y(vh) + γ*(η_ref - η)."""
function continuity_surf_forc!(dη,dUdx,dVdy,η,η_ref,Fη,t)
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

"""Continuity equation's right-hand side with time&space dependent forcing."""
function continuity_forcing!(dη,dUdx,dVdy,η,η_ref,Fη,t)
    m,n = size(dη) .- (2*haloη,2*haloη)
    @boundscheck (m,n+2) == size(dUdx) || throw(BoundsError())
    @boundscheck (m+2,n) == size(dVdy) || throw(BoundsError())
    @boundscheck (m,n) == size(Fη) || throw(BoundsError())

    # avoid recomputation
    Fηtt = Fηt(t)

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            dη[i+1,j+1] = -(dUdx[i,j+1] + dVdy[i+1,j]) + Fηtt*Fη[i,j]
        end
    end
end

"""Continuity equation's right-hand side -∂x(uh) - ∂y(vh) without forcing."""
function continuity_itself!(dη,dUdx,dVdy,η,η_ref,Fη,t)
    m,n = size(dη) .- (2*haloη,2*haloη)
    @boundscheck (m,n+2) == size(dUdx) || throw(BoundsError())
    @boundscheck (m+2,n) == size(dVdy) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            dη[i+1,j+1] = -(dUdx[i,j+1] + dVdy[i+1,j])
        end
    end
end

if surface_relax
    continuity! = continuity_surf_forc!
elseif surface_forcing
    continuity! = continuity_forcing!
else
    continuity! = continuity_itself!
end

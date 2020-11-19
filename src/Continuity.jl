"""Continuity equation's right-hand side with surface relaxation
-∂x(uh) - ∂y(vh) + γ*(η_ref - η)."""
function continuity_surf_relax!(η::Array{T,2},
                                Diag::DiagnosticVars{T,Tprog},
                                S::ModelSetup{T,Tprog},
                                t::Int) where {T,Tprog}

    @unpack dη = Diag.Tendencies
    @unpack dUdx,dVdy = Diag.VolumeFluxes
    @unpack η_ref = S.forcing
    @unpack γ = S.constants

    m,n = size(dη) .- (2,2)
    @boundscheck (m,n+2) == size(dUdx) || throw(BoundsError())
    @boundscheck (m+2,n) == size(dVdy) || throw(BoundsError())
    @boundscheck (m+2,n+2) == size(η) || throw(BoundsError())
    @boundscheck (m,n) == size(η_ref) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            dη[i+1,j+1] = -(Tprog(dUdx[i,j+1]) + Tprog(dVdy[i+1,j])) + Tprog(γ*(η_ref[i,j]-η[i+1,j+1]))
        end
    end
end

"""Continuity equation's right-hand side with time&space dependent forcing."""
function continuity_forcing!(   Diag::DiagnosticVars{T,Tprog},
                                S::ModelSetup{T,Tprog},
                                t::Int) where {T,Tprog}

    @unpack dη = Diag.Tendencies
    @unpack dUdx,dVdy = Diag.VolumeFluxes
    @unpack Fη = S.forcing

    m,n = size(dη) .- (2,2)         # cut off halo
    @boundscheck (m,n+2) == size(dUdx) || throw(BoundsError())
    @boundscheck (m+2,n) == size(dVdy) || throw(BoundsError())
    @boundscheck (m,n) == size(Fη) || throw(BoundsError())

    # avoid recomputation
    @unpack ωFη = S.constants
    Fηt = Ftime(T,t,ωFη)

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            dη[i+1,j+1] = -(Tprog(dUdx[i,j+1]) + Tprog(dVdy[i+1,j])) + Tprog(Fηt*Fη[i,j])
        end
    end
end

"""Continuity equation's right-hand side -∂x(uh) - ∂y(vh) without forcing."""
function continuity_itself!(Diag::DiagnosticVars{T,Tprog},
                            S::ModelSetup{T,Tprog},
                            t::Int) where {T,Tprog}

    @unpack dη = Diag.Tendencies
    @unpack dUdx,dVdy = Diag.VolumeFluxes

    m,n = size(dη) .- (2,2)     # cut off halo
    @boundscheck (m,n+2) == size(dUdx) || throw(BoundsError())
    @boundscheck (m+2,n) == size(dVdy) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            dη[i+1,j+1] = -(Tprog(dUdx[i,j+1]) + Tprog(dVdy[i+1,j]))
        end
    end
end


"""Transit function to call the specified continuity function."""
function continuity!(   u::AbstractMatrix,
                        v::AbstractMatrix,
                        η::AbstractMatrix,
                        Diag::DiagnosticVars,
                        S::ModelSetup,
                        t::Int)
    
    @unpack U,V,dUdx,dVdy = Diag.VolumeFluxes
    @unpack nstep_advcor = S.grid
    @unpack time_scheme,surface_relax,surface_forcing = S.parameters

    if nstep_advcor > 0 ||              # then UV wasn't calculated inside rhs!
        time_scheme != "RK"             # then u,v changed before evaluating semi-implciti continuity
        UVfluxes!(u,v,η,Diag,S)
    end

    # divergence of mass flux
    ∂x!(dUdx,U)
    ∂y!(dVdy,V)

    if surface_relax
        continuity_surf_relax!(η,Diag,S,t)
    elseif surface_forcing
        continuity_forcing!(Diag,S,t)
    else
        continuity_itself!(Diag,S,t)
    end
end
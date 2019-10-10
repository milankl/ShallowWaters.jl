mutable struct Constants{T<:AbstractFloat}

    # RUNGE-KUTTA COEFFICIENTS 3rd/4th order including timestep Δt
    RKaΔt::Array{T,1}
    RKbΔt::Array{T,1}

    # BOUNDARY CONDITIONS
    one_minus_α::T      # tangential boundary condition for the ghost-point copy

    # PHYSICAL CONSTANTS
    g::T                    # gravity
    cD::T                   # quadratic bottom friction - incl grid spacing
    rD::T                   # linear bottom friction - incl grid spacing
    γ::T                    # frequency of interface relaxation
    cSmag::T                # Smagorinsky constant
    νB::T                   # biharmonic diffusion coefficient
    rSST::T                 # tracer restoring timescale
    jSST::T                 # tracer consumption timescale
    SSTmin::T               # tracer minimum

    Constants{T}() where T = new{T}()
end

"""Generator function for the mutable struct Constants."""
function Constants{T}(P::Parameter,G::Grid) where {T<:AbstractFloat}

    C = Constants{T}()

    # Runge-Kutta 3rd/4th order coefficients including time step Δt
    # (which includes the grid spacing Δ too)
    if P.RKo == 3     # version 2
        C.RKaΔt = T.([1/4,0.,3/4]*G.Δt)
        C.RKbΔt = T.([1/3,2/3]*G.Δt)
    elseif P.RKo == 4
        C.RKaΔt = T.([1/6,1/3,1/3,1/6]*G.Δt)
        C.RKbΔt = T.([.5,.5,1.]*G.Δt)
    end

    # for the ghost point copy/tangential boundary conditions
    C.one_minus_α = T(1-P.α)

    C.g = T(P.g)                  # gravity - for Bernoulli potential

    # BOTTOM FRICTION COEFFICENTS
    # incl grid spacing Δ for non-dimensional gradients
    C.cD = T(-G.Δ*P.cD)             # quadratic drag [m]
    C.rD = T(-G.Δ/(P.τD*24*3600))   # linear drag [m/s]

    # INTERFACE RELAXATION FREQUENCY
    # incl grid spacing Δ for non-dimensional gradients
    C.γ = T(G.Δ/(P.t_relax*3600*24))    # [m/s]

    # BIHARMONIC DIFFUSION
    C.cSmag = T(-P.cSmag)   # Smagorinsky coefficient
    C.νB = T(-P.νB/30000)   # linear scaling based on 540m^s/s at Δ=30km

    # TRACER ADVECTION
    C.rSST = T(G.dtadvint/(P.τSST*3600*24))    # tracer restoring [1]
    C.jSST = T(G.dtadvint/(P.jSST*3600*24))    # tracer consumption [1]
    C.SSTmin = T(P.SSTmin)

    return C
end

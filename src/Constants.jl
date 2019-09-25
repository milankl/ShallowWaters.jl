mutable struct Constants{T<:AbstractFloat}

    # RUNGE-KUTTA COEFFICIENTS 3rd/4th order
    RKa::Array{T,1}
    RKb::Array{T,1}

    # BOUNDARY CONDITIONS
    one_minus_α::T      # tangential boundary condition for the ghost-point copy

    # NUMBERS
    zeero::T
    oone::T
    minus_4::T
    one_half::T
    one_twelve::T
    one_quarter::T

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

    # for analysis the old operators with boundary conditions are kept. Constants for these
    # one_over_Δ::T
    # α::T
    # minus_α::T
    # one_minus_α_half::T
    # α_over_Δ::T
    # one_quarter::T
    # minus_3_minus_α::T

    Constants{T}() where T = new{T}()
end

function Constants(::Type{T},P::Parameter,G::Grid) where {T<:AbstractFloat}

    C = Constants{T}()

    # Runge-Kutta 3rd/4th order coefficients
    if P.RKo == 3     # version 2
        C.RKa = T.([1/4,0.,3/4])
        C.RKb = T.([1/3,2/3])
    elseif P.RKo == 4
        C.RKa = T.([1/6,1/3,1/3,1/6])
        C.RKb = T.([.5,.5,1.])
    end

    # for the ghost point copy/tangential boundary conditions
    C.one_minus_α = T(1-P.α)

    C.zeero = zero(T)
    C.oone = one(T)
    C.minus_4 = T(-4.)            # for Laplace operator
    C.one_half = T(0.5)           # for interpolations
    C.one_twelve = T(1/12)
    C.one_quarter = T(0.25)

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

    # for analysis the old operators with boundary conditions are kept. Constants for these
    # C.one_over_Δ = T(1/G.Δ)
    # C.α = T(P.lbc)
    # C.minus_α = T(-P.lbc)
    # C.one_minus_α_half = T(1-0.5*P.lbc)
    # C.α_over_Δ = T(P.lbc/G.Δ)
    # C.one_quarter = T(0.25)
    # C.minus_3_minus_α = T(-3-P.lbc)

    return C
end

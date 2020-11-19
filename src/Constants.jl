struct Constants{T<:AbstractFloat,Tprog<:AbstractFloat}

    # RUNGE-KUTTA COEFFICIENTS 2nd/3rd/4th order including timestep Δt
    RKaΔt::Array{Tprog,1}
    RKbΔt::Array{Tprog,1}
    Δt_Δs::Tprog            # Δt/(s-1) wher s the number of stages
    Δt_Δ::Tprog             # Δt/Δ - timestep divided by grid spacing
    Δt_Δ_half::Tprog        # 1/2 * Δt/Δ

    # BOUNDARY CONDITIONS
    one_minus_α::Tprog      # tangential boundary condition for the ghost-point copy

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
    ωFη::Float64            # frequency [1/s] of seasonal surface forcing incl 2π
    ωFx::Float64            # frequency [1/s] of seasonal wind x incl 2π
    ωFy::Float64            # frequency [1/2] of seasonal wind y incl 2π
end

"""Generator function for the mutable struct Constants."""
function Constants{T,Tprog}(P::Parameter,G::Grid) where {T<:AbstractFloat,Tprog<:AbstractFloat}

    # Runge-Kutta 2nd/3rd/4th order coefficients including time step Δt and grid spacing Δ
    # a are the coefficents to sum the rhs on the fly, such that sum=1
    # b are the coefficents for the update that used for a new evaluation of the RHS
    if P.RKo == 2     # Heun's method
        RKaΔt = Tprog.([1/2,1/2]*G.dtint/G.Δ)
        RKbΔt = Tprog.([1]*G.dtint/G.Δ)
    elseif P.RKo == 3     # version 2 / Heun's 3rd order
        RKaΔt = Tprog.([1/4,0.,3/4]*G.dtint/G.Δ)
        RKbΔt = Tprog.([1/3,2/3]*G.dtint/G.Δ)
    elseif P.RKo == 4
        RKaΔt = Tprog.([1/6,1/3,1/3,1/6]*G.dtint/G.Δ)
        RKbΔt = Tprog.([.5,.5,1.]*G.dtint/G.Δ)
    end

    # Δt/(s-1) for SSPRK2
    Δt_Δs = Tprog(G.dtint/G.Δ/(P.RKs-1))

    # time step and half the time step including the grid spacing as this is not included in the RHS
    Δt_Δ = Tprog(G.dtint/G.Δ)
    Δt_Δ_half = Tprog(G.dtint/G.Δ/2)

    one_minus_α = Tprog(1-P.α)    # for the ghost point copy/tangential boundary conditions
    g = T(P.g)                # gravity - for Bernoulli potential

    # BOTTOM FRICTION COEFFICENTS
    # incl grid spacing Δ for non-dimensional gradients
    cD = T(-G.Δ*P.cD)             # quadratic drag [m]
    rD = T(-G.Δ/(P.τD*24*3600))   # linear drag [m/s]

    # INTERFACE RELAXATION FREQUENCY
    # incl grid spacing Δ for non-dimensional gradients
    γ = T(G.Δ/(P.t_relax*3600*24))    # [m/s]

    # BIHARMONIC DIFFUSION
    cSmag = T(-P.cSmag)   # Smagorinsky coefficient
    νB = T(-P.νB/30000)   # linear scaling based on 540m^s/s at Δ=30km

    # TRACER ADVECTION
    rSST = T(G.dtadvint/(P.τSST*3600*24))    # tracer restoring [1]
    jSST = T(G.dtadvint/(P.jSST*3600*24))    # tracer consumption [1]
    SSTmin = T(P.SSTmin)

    # TIME DEPENDENT FORCING
    ωFη = -2π*P.ωFη/24/365.25/3600
    ωFx = 2π*P.ωFx/24/365.25/3600
    ωFy = 2π*P.ωFy/24/365.25/3600

    return Constants{T,Tprog}(  RKaΔt,RKbΔt,Δt_Δs,Δt_Δ,Δt_Δ_half,
                                one_minus_α,
                                g,cD,rD,γ,cSmag,νB,rSST,
                                jSST,SSTmin,ωFη,ωFx,ωFy)
end

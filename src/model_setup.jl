mutable struct PrognosticVars{T<:AbstractFloat}
    u::Array{T,2}           # u-velocity
    v::Array{T,2}           # v-velocity
    η::Array{T,2}           # sea surface height / interface displacement
    sst::Array{T,2}         # tracer / sea surface temperature
end

struct DiagnosticVars{T,Tprog}
    RungeKutta::RungeKuttaVars{Tprog}
    Tendencies::TendencyVars{Tprog}
    VolumeFluxes::VolumeFluxVars{T}
    Vorticity::VorticityVars{T}
    Bernoulli::BernoulliVars{T}
    Bottomdrag::BottomdragVars{T}
    ArakawaHsu::ArakawaHsuVars{T}
    Laplace::LaplaceVars{T}
    Smagorinsky::SmagorinskyVars{T}
    SemiLagrange::SemiLagrangeVars{T}
    PrognosticVarsRHS::PrognosticVars{T}        # low precision version
end

mutable struct ModelSetup{T<:AbstractFloat,Tprog<:AbstractFloat}
    parameters::Parameter
    grid::Grid{T,Tprog}
    constants::Constants{T,Tprog}
    forcing::Forcing{T}
    Prog::PrognosticVars{Tprog}
    Diag::DiagnosticVars{T, Tprog}
    t::Int
end
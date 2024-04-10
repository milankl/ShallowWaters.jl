"""Runge Kutta time stepping scheme diagnostic cariables collected in a struct."""
@with_kw struct RungeKuttaVars{T<:AbstractFloat}

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = ny-1                                     # v-grid in y-direction
    nqx::Int = if (bc == "periodic") nx else nx+1 end   # q-grid in x-direction
    nqy::Int = ny+1                                     # q-grid in y-direction

    # EDGE POINT (1 = yes, 0 = no)
    ep::Int = if bc == "periodic" 1 else 0 end  # is there a u-point on the left edge?

    u0::Array{T,2} = zeros(T,nux+2*halo,nuy+2*halo)     # u-velocities for RK updates
    u1::Array{T,2} = zeros(T,nux+2*halo,nuy+2*halo)
    v0::Array{T,2} = zeros(T,nvx+2*halo,nvy+2*halo)     # v-velocities for RK updates
    v1::Array{T,2} = zeros(T,nvx+2*halo,nvy+2*halo)
    η0::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη)     # sea surface height for RK updates
    η1::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη)
end

"""Generator function for RungeKutta VarCollection."""
function RungeKuttaVars{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return RungeKuttaVars{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

###################################################

"""Tendencies collected in a struct."""
@with_kw struct TendencyVars{T<:AbstractFloat}

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = ny-1                                     # v-grid in y-direction
    nqx::Int = if (bc == "periodic") nx else nx+1 end   # q-grid in x-direction
    nqy::Int = ny+1                                     # q-grid in y-direction

    # EDGE POINT (1 = yes, 0 = no)
    ep::Int = if bc == "periodic" 1 else 0 end  # is there a u-point on the left edge?

    du::Array{T,2} = zeros(T,nux+2*halo,nuy+2*halo)     # tendency of u without time step
    dv::Array{T,2} = zeros(T,nvx+2*halo,nvy+2*halo)     # tendency of v without time step
    dη::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη)     # tendency of η without time step

    # sum of tendencies (incl time step) over all sub-steps
    du_sum::Array{T,2} = zeros(T,nux+2*halo,nuy+2*halo) 
    dv_sum::Array{T,2} = zeros(T,nvx+2*halo,nvy+2*halo)
    dη_sum::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη)

    # compensation for tendencies (variant of Kahan summation)
    du_comp::Array{T,2} = zeros(T,nux+2*halo,nuy+2*halo) 
    dv_comp::Array{T,2} = zeros(T,nvx+2*halo,nvy+2*halo)
    dη_comp::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη)
end

"""Generator function for Tendencies VarCollection."""
function TendencyVars{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return TendencyVars{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

###########################################################

"""VolumeFluxes collected in a struct."""
@with_kw struct VolumeFluxVars{T<:AbstractFloat}

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = ny-1                                     # v-grid in y-direction
    nqx::Int = if (bc == "periodic") nx else nx+1 end   # q-grid in x-direction
    nqy::Int = ny+1                                     # q-grid in y-direction

    # EDGE POINT (1 = yes, 0 = no)
    ep::Int = if bc == "periodic" 1 else 0 end      # is there a u-point on the left edge?

    h::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη)         # layer thickness
    h_u::Array{T,2} = zeros(T,nx+2*haloη-1,ny+2*haloη)     # layer thickness on u-grid
    U::Array{T,2} = zeros(T,nx+2*haloη-1,ny+2*haloη)       # U=uh volume flux

    h_v::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη-1)     # layer thickness on v-grid
    V::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη-1)       # V=vh volume flux

    dUdx::Array{T,2} = zeros(T,nx+2*haloη-2,ny+2*haloη)    # gradients thereof
    dVdy::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη-2)
end

"""Generator function for VolumeFluxes VarCollection."""
function VolumeFluxVars{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return VolumeFluxVars{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

###############################################################

"""Vorticity variables collected in a struct."""
@with_kw struct VorticityVars{T<:AbstractFloat}

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = ny-1                                     # v-grid in y-direction
    nqx::Int = if (bc == "periodic") nx else nx+1 end   # q-grid in x-direction
    nqy::Int = ny+1                                     # q-grid in y-direction

    # EDGE POINT (1 = yes, 0 = no)
    ep::Int = if bc == "periodic" 1 else 0 end      # is there a u-point on the left edge?

    h_q::Array{T,2} = zeros(T,nx+2*haloη-1,ny+2*haloη-1)  # layer thickness h interpolated on q-grid
    q::Array{T,2} = zeros(T,nx+2*haloη-1,ny+2*haloη-1)    # potential vorticity

    q_v::Array{T,2} = zeros(T,nx+2*haloη-2,ny+2*haloη-1)  # q interpolated on v-grid
    U_v::Array{T,2} = zeros(T,nx+2*haloη-2,ny+2*haloη-1)  # mass flux U=uh on v-grid

    q_u::Array{T,2} = zeros(T,nx+2*haloη-1,ny+2*haloη-2)  # q interpolated on u-grid
    V_u::Array{T,2} = zeros(T,nx+2*haloη-1,ny+2*haloη-2)  # mass flux V=vh on v-grid

    qhu::Array{T,2} = zeros(T,nvx,nvy)            # potential vorticity advection term u-component
    qhv::Array{T,2} = zeros(T,nux,nuy)            # potential vorticity advection term v-component

    u_v::Array{T,2} = zeros(T,nux+2*halo-1,nuy+2*halo-1)  # u-velocity on v-grid
    v_u::Array{T,2} = zeros(T,nvx+2*halo-1,nvy+2*halo-1)  # v-velocity on u-grid

    dudx::Array{T,2} = zeros(T,nux+2*halo-1,nuy+2*halo)   # ∂u/∂x
    dudy::Array{T,2} = zeros(T,nux+2*halo,nuy+2*halo-1)   # ∂u/∂y

    dvdx::Array{T,2} = zeros(T,nvx+2*halo-1,nvy+2*halo)   # ∂v/∂x
    dvdy::Array{T,2} = zeros(T,nvx+2*halo,nvy+2*halo-1)   # ∂v/∂y
end

"""Generator function for Vorticity VarCollection."""
function VorticityVars{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return VorticityVars{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

####################################################################

"""Bernoulli variables collected in a struct."""
@with_kw struct BernoulliVars{T<:AbstractFloat}

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = ny-1                                     # v-grid in y-direction
    nqx::Int = if (bc == "periodic") nx else nx+1 end   # q-grid in x-direction
    nqy::Int = ny+1                                     # q-grid in y-direction

    # EDGE POINT (1 = yes, 0 = no)
    ep::Int = if bc == "periodic" 1 else 0 end      # is there a u-point on the left edge?

    u²::Array{T,2} = zeros(T,nux+2*halo,nuy+2*halo)         # u-velocity squared
    v²::Array{T,2} = zeros(T,nvx+2*halo,nvy+2*halo)         # v-velocity squared

    KEu::Array{T,2} = zeros(T,nux+2*halo-1,nuy+2*halo)      # u-velocity squared on T-grid
    KEv::Array{T,2} = zeros(T,nvx+2*halo,nvy+2*halo-1)      # v-velocity squared on T-grid

    p::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη)          # Bernoulli potential
    dpdx::Array{T,2} = zeros(T,nx+2*haloη-1,ny+2*haloη)     # ∂p/∂x
    dpdy::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη-1)     # ∂p/∂y
end

"""Generator function for Bernoulli VarCollection."""
function BernoulliVars{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return BernoulliVars{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

####################################################################

"""Bottomdrag variables collected in a struct."""
@with_kw struct BottomdragVars{T<:AbstractFloat}

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = ny-1                                     # v-grid in y-direction
    nqx::Int = if (bc == "periodic") nx else nx+1 end   # q-grid in x-direction
    nqy::Int = ny+1                                     # q-grid in y-direction

    # EDGE POINT (1 = yes, 0 = no)
    ep::Int = if bc == "periodic" 1 else 0 end      # is there a u-point on the left edge?

    sqrtKE::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη)       # sqrt of kinetic energy
    sqrtKE_u::Array{T,2} = zeros(T,nx+2*haloη-1,ny+2*haloη)   # interpolated on u-grid
    sqrtKE_v::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη-1)   # interpolated on v-grid

    Bu::Array{T,2} = zeros(T,nx+2*haloη-1,ny+2*haloη)         # bottom friction term u-component
    Bv::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη-1)         # bottom friction term v-component
end

"""Generator function for Bottomdrag VarCollection."""
function BottomdragVars{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return BottomdragVars{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

####################################################################

"""ArakawaHsu variables collected in a struct."""
@with_kw struct ArakawaHsuVars{T<:AbstractFloat}

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = ny-1                                     # v-grid in y-direction
    nqx::Int = if (bc == "periodic") nx else nx+1 end   # q-grid in x-direction
    nqy::Int = ny+1                                     # q-grid in y-direction

    # EDGE POINT (1 = yes, 0 = no)
    ep::Int = if bc == "periodic" 1 else 0 end      # is there a u-point on the left edge?

    # Linear combination of potential vorticity
    qα::Array{T,2} = zeros(T,nx+2*haloη-2,ny+2*haloη-2)
    qβ::Array{T,2} = zeros(T,nx+2*haloη-1,ny+2*haloη-2)
    qγ::Array{T,2} = zeros(T,nx+2*haloη-1,ny+2*haloη-2)
    qδ::Array{T,2} = zeros(T,nx+2*haloη-2,ny+2*haloη-2)
end

"""Generator function for ArakawaHsu VarCollection."""
function ArakawaHsuVars{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return ArakawaHsuVars{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

####################################################################

"""Laplace variables collected in a struct."""
@with_kw struct LaplaceVars{T<:AbstractFloat}

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = ny-1                                     # v-grid in y-direction
    nqx::Int = if (bc == "periodic") nx else nx+1 end   # q-grid in x-direction
    nqy::Int = ny+1                                     # q-grid in y-direction

    # EDGE POINT (1 = yes, 0 = no)
    ep::Int = if bc == "periodic" 1 else 0 end      # is there a u-point on the left edge?

    Lu::Array{T,2} = zeros(T,nux+2*halo-2,nuy+2*halo-2)         # ∇²u
    Lv::Array{T,2} = zeros(T,nvx+2*halo-2,nvy+2*halo-2)         # ∇²v

    # Derivatives of Lu,Lv
    dLudx::Array{T,2} = zeros(T,nux+2*halo-3,nuy+2*halo-2)
    dLudy::Array{T,2} = zeros(T,nux+2*halo-2,nuy+2*halo-3)
    dLvdx::Array{T,2} = zeros(T,nvx+2*halo-3,nvy+2*halo-2)
    dLvdy::Array{T,2} = zeros(T,nvx+2*halo-2,nvy+2*halo-3)
end

"""Generator function for Laplace VarCollection."""
function LaplaceVars{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return LaplaceVars{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

####################################################################

"""Smagorinsky variables collected in a struct."""
@with_kw struct SmagorinskyVars{T<:AbstractFloat}

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = ny-1                                     # v-grid in y-direction
    nqx::Int = if (bc == "periodic") nx else nx+1 end   # q-grid in x-direction
    nqy::Int = ny+1                                     # q-grid in y-direction

    # EDGE POINT (1 = yes, 0 = no)
    ep::Int = if bc == "periodic" 1 else 0 end      # is there a u-point on the left edge?

    DT::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη)       # Tension squared (on the T-grid)
    DS::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη)       # Shearing strain squared (on the T-grid)
    νSmag::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη)    # Viscosity coefficient

    # Tension squared on the q-grid
    DS_q::Array{T,2} = zeros(T,nvx+2*halo-1,nvy+2*halo)

    # Smagorinsky viscosity coefficient on the q-grid
    νSmag_q::Array{T,2} = zeros(T,nx+2*haloη-1,ny+2*haloη-1)

    # Entries of the Smagorinsky viscous tensor
    S12::Array{T,2} = zeros(T,nx+2*haloη-1,ny+2*haloη-1)
    S21::Array{T,2} = zeros(T,nx+2*haloη-1,ny+2*haloη-1)

    S11::Array{T,2} = zeros(T,nux+2*halo-3,nuy+2*halo-2)
    S22::Array{T,2} = zeros(T,nvx+2*halo-2,nvy+2*halo-3)

    # u- and v-components 1 and 2 of the biharmonic diffusion tendencies
    LLu1::Array{T,2} = zeros(T,nux+2*halo-4,nuy+2*halo-2)
    LLu2::Array{T,2} = zeros(T,nx+1,ny)

    LLv1::Array{T,2} = zeros(T,nx,ny+1)
    LLv2::Array{T,2} = zeros(T,nvx+2*halo-2,nvy+2*halo-4)
end

"""Generator function for Smagorinsky VarCollection."""
function SmagorinskyVars{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return SmagorinskyVars{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

####################################################################

"""SemiLagrange variables collected in a struct."""
@with_kw struct SemiLagrangeVars{T<:AbstractFloat}

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int
    halosstx::Int
    halossty::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = ny-1                                     # v-grid in y-direction
    nqx::Int = if (bc == "periodic") nx else nx+1 end   # q-grid in x-direction
    nqy::Int = ny+1                                     # q-grid in y-direction

    # EDGE POINT (1 = yes, 0 = no)
    ep::Int = if bc == "periodic" 1 else 0 end      # is there a u-point on the left edge?

    xd::Array{T,2} = zeros(T,nx,ny)                         # departure points x-coord
    yd::Array{T,2} = zeros(T,nx,ny)                         # departure points y-coord

    um::Array{T,2} = zeros(T,nux+2*halo,nuy+2*halo)         # u-velocity temporal mid-point
    vm::Array{T,2} = zeros(T,nvx+2*halo,nvy+2*halo)         # v-velocity temporal mid-point

    u_T::Array{T,2} = zeros(T,nux+2*halo-1,nuy+2*halo)      # u-velocity interpolated on T-grid
    um_T::Array{T,2} = zeros(T,nux+2*halo-1,nuy+2*halo)     # um interpolated on T-grid
    v_T::Array{T,2} = zeros(T,nvx+2*halo,nvy+2*halo-1)      # v-velocity interpolated on T-grid
    vm_T::Array{T,2} = zeros(T,nvx+2*halo,nvy+2*halo-1)     # vm interpolated on T-grid

    uinterp::Array{T,2} = zeros(T,nx,ny)                    # u interpolated on mid-point xd,yd
    vinterp::Array{T,2} = zeros(T,nx,ny)                    # v interpolated on mid-point xd,yd

    ssti::Array{T,2} = zeros(T,nx+2*halosstx,ny+2*halossty) # sst interpolated on departure points
    sst_ref::Array{T,2} = zeros(T,nx+2*halosstx,ny+2*halossty) # sst initial conditions for relaxation

    # compensated summation
    dsst_comp::Array{T,2} = zeros(T,nx+2*halosstx,ny+2*halossty)
end

"""Generator function for SemiLagrange VarCollection."""
function SemiLagrangeVars{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G
    @unpack halosstx,halossty = G

    return SemiLagrangeVars{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη,
                            halosstx=halosstx,halossty=halossty)
end

###################################################################

# struct DiagnosticVars{T,Tprog}
#     RungeKutta::RungeKuttaVars{Tprog}
#     Tendencies::TendencyVars{Tprog}
#     VolumeFluxes::VolumeFluxVars{T}
#     Vorticity::VorticityVars{T}
#     Bernoulli::BernoulliVars{T}
#     Bottomdrag::BottomdragVars{T}
#     ArakawaHsu::ArakawaHsuVars{T}
#     Laplace::LaplaceVars{T}
#     Smagorinsky::SmagorinskyVars{T}
#     SemiLagrange::SemiLagrangeVars{T}
#     PrognosticVarsRHS::PrognosticVars{T}        # low precision version
# end

"""Preallocate the diagnostic variables and return them as matrices in structs."""
function preallocate(   ::Type{T},
                        ::Type{Tprog},
                        G::Grid) where {T<:AbstractFloat,Tprog<:AbstractFloat}

    RK = RungeKuttaVars{Tprog}(G)
    TD = TendencyVars{Tprog}(G)
    VF = VolumeFluxVars{T}(G)
    VT = VorticityVars{T}(G)
    BN = BernoulliVars{T}(G)
    BD = BottomdragVars{T}(G)
    AH = ArakawaHsuVars{T}(G)
    LP = LaplaceVars{T}(G)
    SM = SmagorinskyVars{T}(G)
    SL = SemiLagrangeVars{T}(G)
    PV = PrognosticVars{T}(G)

    return DiagnosticVars{T,Tprog}(RK,TD,VF,VT,BN,BD,AH,LP,SM,SL,PV)
end

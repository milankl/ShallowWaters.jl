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

""" Variables that appear in Zanna-Bolton forcing term """
@with_kw struct ZBVars{T<:AbstractFloat}

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int
    halosstx::Int
    halossty::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end      # u-grid in x-direction
    nuy::Int = ny                                          # u-grid in y-direction
    nvx::Int = nx                                          # v-grid in x-direction
    nvy::Int = ny-1                                        # v-grid in y-direction
    nqx::Int = if (bc == "periodic") nx else nx+1 end      # q-grid in x-direction
    nqy::Int = ny+1                                        # q-grid in y-direction

    dudx::Array{T,2} = zeros(T,nux+2*halo-1,nuy+2*halo)    # ∂u/∂x
    dudy::Array{T,2} = zeros(T,nux+2*halo,nuy+2*halo-1)    # ∂u/∂y
    dvdx::Array{T,2} = zeros(T,nvx+2*halo-1,nvy+2*halo)    # ∂v/∂x
    dvdy::Array{T,2} = zeros(T,nvx+2*halo,nvy+2*halo-1)    # ∂v/∂y

    γ₀::T=0.3                       # coefficient in parameterization term

    ζ::Array{T,2} = zeros(T,nqx,nqy)      # relative vorticity, cell corners
    ζsq::Array{T,2} = zeros(T,nqx,nqy)    # relative vorticity squared, cell corners

    D::Array{T,2} = zeros(T,nqx,nqy)      # shear deformation of flow field, cell corners
    Dsq::Array{T,2} = zeros(T,nqx,nqy)    # square of the tensor

    D_n::Array{T,2} = zeros(T,nvx+2*halo-1,nvy+2*halo)
    D_nT::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη) 
    D_q::Array{T,2} = zeros(T,nqx,nqy)

    Dhat::Array{T,2} = zeros(T,nqx-1+2*haloη,nqy-1+2*haloη)     # stretch deformation of flow field, cell centers w/ halo
    Dhatsq::Array{T,2} = zeros(T,nqx-1+2*haloη,nqy-1+2*haloη)   # square of the tensor
    Dhatq::Array{T,2} = zeros(T,nqx,nqy)                  # tensor interpolated onto q-grid

    ζsqT::Array{T,2} = zeros(T,nqx-1,nqy-1)     # ζ^2 interpolated to cell centers
    ζD::Array{T,2} = zeros(T,nqx,nqy)           # ζD, cell corners
    ζDT::Array{T,2} = zeros(T,nqx-1,nqy-1)      # ζD, placed on cell centers
    ζDhat::Array{T,2} = zeros(T,nqx,nqy)        # ζDhat, cell corners
    
    trace::Array{T,2} = zeros(T,nx,ny)     # ζ^2 (+ D^2 + Dhat^2), cell centers. We only compute ζ^2 rather than the whole sum

    ζD_filtered::Array{T,2} = zeros(T,nqx-1,nqy-1)        # ζD with filter applied
    ζDhat_filtered::Array{T,2} = zeros(T,nqx,nqy)         # ζDhat with filter applied
    trace_filtered::Array{T,2} = zeros(T,nqx-1,nqy-1)     # trace with filter applied

    dζDdx::Array{T,2} = zeros(T,nux,nuy)             # u-grid
    dζDhatdy::Array{T,2} = zeros(T,nux+halo,nuy)     # u-grid, initially with extra halo points
    dtracedx::Array{T,2} = zeros(T,nux,nuy)          # u-grid 

    S_u::Array{T,2} = zeros(T,nux,nuy)             # total forcing in x-direction

    dζDhatdx::Array{T,2} = zeros(T,nvx,nvy+halo)   # v-grid, initially with extra halo points
    dζDdy::Array{T,2} = zeros(T,nvx,nvy)           # v-grid
    dtracedy::Array{T,2} = zeros(T,nvx,nvy)        # v-grid

    S_v::Array{T,2} = zeros(T,nvx,nvy)             # total forcing in y-direction

end

"""Generator function for ZB_momentum terms."""
function ZBVars{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G
    @unpack halosstx,halossty = G

    return ZBVars{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη,
                            halosstx=halosstx,halossty=halossty)
end

""" Variables that appear in NN forcing term """
@with_kw mutable struct CNNVars{T<:AbstractFloat, SuLayerType, SvLayerType, SuModelType, SvModelType}#, SuCompiledType, SvCompiledType, DSuCompiledType, DSvCompiledType}

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int
    halosstx::Int
    halossty::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end      # u-grid in x-direction
    nuy::Int = ny                                          # u-grid in y-direction
    nvx::Int = nx                                          # v-grid in x-direction
    nvy::Int = ny-1                                        # v-grid in y-direction
    nqx::Int = if (bc == "periodic") nx else nx+1 end      # q-grid in x-direction
    nqy::Int = ny+1                                        # q-grid in y-direction

    γ₀::Float64=0.3                       # coefficient in parameterization term

    dudx::Array{T,2} = zeros(T,nux+2*halo-1,nuy+2*halo)    # ∂u/∂x
    dudy::Array{T,2} = zeros(T,nux+2*halo,nuy+2*halo-1)    # ∂u/∂y
    dvdx::Array{T,2} = zeros(T,nvx+2*halo-1,nvy+2*halo)    # ∂v/∂x
    dvdy::Array{T,2} = zeros(T,nvx+2*halo,nvy+2*halo-1)    # ∂v/∂y

    ζ::Array{T,2} = zeros(T,nqx,nqy)      # relative vorticity, cell corners
    D::Array{T,2} = zeros(T,nqx,nqy)      # shear deformation of flow field, cell corners
    Dhat::Array{T,2} = zeros(T,nqx-1+2*haloη,nqy-1+2*haloη)     # stretch deformation of flow field, cell centers w/ halo

    Dhatq::Array{T,2} = zeros(T,nqx,nqy)    # stretch deformation, interpolated to cell corners to match ζ and D

    ζT::Array{T,2} = zeros(T,nqx-1,nqy-1)         # ζ interpolated to cell centers
    DT::Array{T,2} = zeros(T,nqx-1,nqy-1)         # D, interpolated on cell centers
    DhatT::Array{T,2} = zeros(T,nqx-1,nqy-1)      # Dhat, further interpolated to cell centers, now with no halo

    T11::Array{T,2} = zeros(T,nx,ny)
    T12::Array{T,2} = zeros(T,nqx,nqy)
    T22::Array{T,2} = zeros(T,nx,ny)

    dT11dx::Array{T,2} = zeros(T,nux,nuy)    # derivative of T11 in the x-direction, u-grid
    dT12dy::Array{T,2} = zeros(T,nux+halo,nuy)    # derivative of T12 in the y-direction, u-grid
    dT12dx::Array{T,2} = zeros(T,nvx,nvy+halo)    # derivative of T12 in the x-direction, v-grid
    dT22dy::Array{T,2} = zeros(T,nvx,nvy)    # derivative of T22 in the y-direction, v-grid

    res_Su::Array{T,2} = zeros(nqx,nuy)
    res_Sv::Array{T,2} = zeros(nvx,nqy)

    S_u::Array{T,2} = zeros(T,nux,nuy)             # total forcing in x-direction
    S_v::Array{T,2} = zeros(T,nvx,nvy)             # total forcing in y-direction

    Su_layers::SuLayerType
    Sv_layers::SvLayerType

    model_Su::SuModelType
    model_Sv::SvModelType

    # compiled_Su::SuCompiledType
    # compiled_Sv::SvCompiledType

    # compiled_dSu::DSuCompiledType
    # compiled_dSv::DSvCompiledType

end

"""Generator function for convolutional NN momentum terms"""
function CNNVars{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc,Δ= G
    @unpack halo,haloη = G
    @unpack halosstx,halossty = G

    nqx = if (bc == "periodic") nx else nx+1 end      # q-grid in x-direction
    nqy = ny+1                                        # q-grid in y-direction

    # This was the size of the CNNs set for my work. There's currently no setup for the user
    # to decide how large/small to make the CNN forcing term, the only way to alter the number of
    # weights is to manually change these values
    Su_dims = [3,25,25,1]
    Sv_dims = [3,25,25,2]

    Su_layers = Lux.Chain(
        (
            Lux.Conv((5,5), Su_dims[i] => Su_dims[i+1], (i == (length(Su_dims)-1) ? identity : gelu); pad=SamePad(),use_bias=false)
            for i in 1:(length(Su_dims)-1)
        )...
    )

    Sv_layers = Lux.Chain(
        (
            Lux.Conv((5,5), Sv_dims[i] => Sv_dims[i+1], (i == (length(Sv_dims)-1) ? identity : gelu); pad=SamePad(),use_bias=false)
            for i in 1:(length(Sv_dims)-1)
        )...
    )

    model_Su = Lux.setup(Random.default_rng(), Su_layers)
    model_Sv = Lux.setup(Random.default_rng(), Sv_layers)

    use_reactant = false
    # if use_reactant
    #     model_Su = Reactant.to_rarray(model_Su)
    #     Su_input = Reactant.to_rarray(Array{T}(undef, 9+9+4, nqx, nqy))
    #     Sv_input = Reactant.to_rarray(Array{T}(undef, 9+4+4, nx, ny))

    #     Su_dinput = Reactant.to_rarray(Array{T}(undef, 9+9+4, nqx, nqy))
    #     Sv_dinput = Reactant.to_rarray(Array{T}(undef, 9+4+4, nx, ny))

    #     d_Su_res = Reactant.to_rarray(Array{T}(undef, 1, nqx, nqy))
    #     d_Sv_res = Reactant.to_rarray(Array{T}(undef, 2, nx, ny))
    # end
    # if use_reactant
    #     model_Sv = Reactant.to_rarray(model_Sv)
    # end

    # if use_reactant
    #     compiled_Su = Reactant.@compile Lux.apply(Su_layers, Su_input, model_Su[1], model_Su[2])
    #     compiled_Sv = Reactant.@compile Lux.apply(Sv_layers, Sv_input, model_Sv[1], model_Sv[2])

    #     compiled_dSu = Reactant.@compile grad_apply(d_Su_res, deepcopy(model_Su[1]), Su_layers, Su_input, Su_dinput, model_Su[1], model_Su[2])
    #     compiled_dSv = Reactant.@compile grad_apply(d_Sv_res, deepcopy(model_Sv[1]), Sv_layers, Sv_input, Sv_dinput, model_Sv[1], model_Sv[2])
    # else
    #     compiled_Su = nothing
    #     compiled_Sv = nothing
    #     compiled_dSu = nothing
    #     compiled_dSv = nothing
    # end

    return CNNVars{T, typeof(Su_layers), typeof(Sv_layers), typeof(model_Su), typeof(model_Sv)}(; nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη,
                    halosstx=halosstx,halossty=halossty, Su_layers, Sv_layers, model_Su, model_Sv#, compiled_Su, compiled_Sv, compiled_dSu, compiled_dSv
    )
end

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
    ZB = ZBVars{Tprog}(G)
    CNN = CNNVars{Tprog}(G)

    return DiagnosticVars{T,Tprog}(RK,TD,VF,VT,BN,BD,AH,LP,SM,SL,PV,ZB,CNN)
end

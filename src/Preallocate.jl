abstract type VarCollection end

struct DiagnosticVars
    RungeKutta::VarCollection
    Tendencies::VarCollection
    VolumeFluxes::VarCollection
    Vorticity::VarCollection
    Bernoulli::VarCollection
    Bottomdrag::VarCollection
    ArakawaHsu::VarCollection
    Laplace::VarCollection
    Smagorinsky::VarCollection
    SemiLagrange::VarCollection
end

"""Preallocate the diagnostic variables and return them as matrices in structs."""
function preallocate(::Type{T},G::Grid) where {T<:AbstractFloat}
    RK = RungeKutta{T}(G)
    TD = Tendencies{T}(G)
    VF = VolumeFluxes{T}(G)
    VT = Vorticity{T}(G)
    BN = Bernoulli{T}(G)
    BD = Bottomdrag{T}(G)
    AH = ArakawaHsu{T}(G)
    LP = Laplace{T}(G)
    SM = Smagorinsky{T}(G)
    SL = SemiLagrange{T}(G)

    return DiagnosticVars(RK,TD,VF,VT,BN,BD,AH,LP,SM,SL)
end

"""Runge Kutta time stepping scheme diagnostic cariables collected in a struct."""
@with_kw struct RungeKutta{T<:AbstractFloat} <: VarCollection

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = nx-1                                     # v-grid in y-direction
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
function RungeKutta{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return RungeKutta{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

###################################################

"""Tendencies collected in a struct."""
@with_kw struct Tendencies{T<:AbstractFloat} <: VarCollection

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = nx-1                                     # v-grid in y-direction
    nqx::Int = if (bc == "periodic") nx else nx+1 end   # q-grid in x-direction
    nqy::Int = ny+1                                     # q-grid in y-direction

    # EDGE POINT (1 = yes, 0 = no)
    ep::Int = if bc == "periodic" 1 else 0 end  # is there a u-point on the left edge?

    du::Array{T,2} = zeros(T,nux+2*halo,nuy+2*halo)     # tendency of u without time step
    dv::Array{T,2} = zeros(T,nvx+2*halo,nvy+2*halo)     # tendency of v without time step
    dη::Array{T,2} = zeros(T,nx+2*haloη,ny+2*haloη)     # tendency of η without time step
end

"""Generator function for Tendencies VarCollection."""
function Tendencies{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return Tendencies{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

###########################################################

"""VolumeFluxes collected in a struct."""
@with_kw struct VolumeFluxes{T<:AbstractFloat} <: VarCollection

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = nx-1                                     # v-grid in y-direction
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
function VolumeFluxes{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return VolumeFluxes{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

###############################################################

"""Vorticity variables collected in a struct."""
@with_kw struct Vorticity{T<:AbstractFloat} <: VarCollection

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = nx-1                                     # v-grid in y-direction
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
function Vorticity{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return Vorticity{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

####################################################################

"""Bernoulli variables collected in a struct."""
@with_kw struct Bernoulli{T<:AbstractFloat} <: VarCollection

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = nx-1                                     # v-grid in y-direction
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
function Bernoulli{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return Bernoulli{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

####################################################################

"""Bottomdrag variables collected in a struct."""
@with_kw struct Bottomdrag{T<:AbstractFloat} <: VarCollection

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = nx-1                                     # v-grid in y-direction
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
function Bottomdrag{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return Bottomdrag{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

####################################################################

"""ArakawaHsu variables collected in a struct."""
@with_kw struct ArakawaHsu{T<:AbstractFloat} <: VarCollection

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = nx-1                                     # v-grid in y-direction
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
function ArakawaHsu{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return ArakawaHsu{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

####################################################################

"""Laplace variables collected in a struct."""
@with_kw struct Laplace{T<:AbstractFloat} <: VarCollection

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = nx-1                                     # v-grid in y-direction
    nqx::Int = if (bc == "periodic") nx else nx+1 end   # q-grid in x-direction
    nqy::Int = ny+1                                     # q-grid in y-direction

    # EDGE POINT (1 = yes, 0 = no)
    ep::Int = if bc == "periodic" 1 else 0 end      # is there a u-point on the left edge?

    Lu = zeros(T,nux+2*halo-2,nuy+2*halo-2)         # ∇²u
    Lv = zeros(T,nvx+2*halo-2,nvy+2*halo-2)         # ∇²v

    # Derivatives of Lu,Lv
    dLudx = zeros(T,nux+2*halo-3,nuy+2*halo-2)
    dLudy = zeros(T,nux+2*halo-2,nuy+2*halo-3)
    dLvdx = zeros(T,nvx+2*halo-3,nvy+2*halo-2)
    dLvdy = zeros(T,nvx+2*halo-2,nvy+2*halo-3)
end

"""Generator function for Laplace VarCollection."""
function Laplace{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return Laplace{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

####################################################################

"""Smagorinsky variables collected in a struct."""
@with_kw struct Smagorinsky{T<:AbstractFloat} <: VarCollection

    # to be specified
    nx::Int
    ny::Int
    bc::String
    halo::Int
    haloη::Int

    nux::Int = if (bc == "periodic") nx else nx-1 end   # u-grid in x-direction
    nuy::Int = ny                                       # u-grid in y-direction
    nvx::Int = nx                                       # v-grid in x-direction
    nvy::Int = nx-1                                     # v-grid in y-direction
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
function Smagorinsky{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G

    return Smagorinsky{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη)
end

####################################################################

"""SemiLagrange variables collected in a struct."""
@with_kw struct SemiLagrange{T<:AbstractFloat} <: VarCollection

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
    nvy::Int = nx-1                                     # v-grid in y-direction
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
end

"""Generator function for SemiLagrange VarCollection."""
function SemiLagrange{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,bc = G
    @unpack halo,haloη = G
    @unpack halosstx,halossty = G

    return SemiLagrange{T}(nx=nx,ny=ny,bc=bc,halo=halo,haloη=haloη,
                            halosstx=halosstx,halossty=halossty)
end

abstract type NamedTupleInStruct end

""" Create structs for different sets of diagnostic variables, each containing a NamedTuple with all matrices."""
for XVars in (  :RungeKutta,
                :Tendencies,
                :VolumeFluxes,
                :Vorticity,
                :Bernoulli,
                :Bottomdrag,
                :Sadourny,
                :ArakawaHsu,
                :Laplace,
                :Smagorinsky,
                :SemiLagrange)
    @eval begin
        struct $XVars{NamedTup} <: NamedTupleInStruct
            data::NamedTup
        end
    end
end

"""Propagate the struct.field notation for NamedTupleInStruct."""
Base.getproperty(S::NamedTupleInStruct,field::Symbol) = getfield(getfield(S,:data), field)

struct DiagnosticVars
    RungeKutta::NamedTupleInStruct
    Tendencies::NamedTupleInStruct
    VolumeFluxes::NamedTupleInStruct
    Vorticity::NamedTupleInStruct
    Bernoulli::NamedTupleInStruct
    Bottomdrag::NamedTupleInStruct
    ArakawaHsu::NamedTupleInStruct
    Laplace::NamedTupleInStruct
    Smagorinsky::NamedTupleInStruct
    SemiLagrange::NamedTupleInStruct
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

function RungeKutta{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo,haloη = G

    NT = (  u0 = Array{T,2}(undef,nux+2*halo,nuy+2*halo),       # u-velocities for RK updates
            u1 = Array{T,2}(undef,nux+2*halo,nuy+2*halo),

            v0 = Array{T,2}(undef,nvx+2*halo,nvy+2*halo),       # v-velocities for RK updates
            v1 = Array{T,2}(undef,nvx+2*halo,nvy+2*halo),

            η0 = Array{T,2}(undef,nx+2*haloη,ny+2*haloη),       # sea surface height for RK updates
            η1 = Array{T,2}(undef,nx+2*haloη,ny+2*haloη)
            )
    return RungeKutta{typeof(NT)}(NT)
end

function Tendencies{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo,haloη = G

    NT = (  du = Array{T,2}(undef,nux+2*halo,nuy+2*halo),       # tendency of u without time step
            dv = Array{T,2}(undef,nvx+2*halo,nvy+2*halo),       # tendency of v without time step
            dη = Array{T,2}(undef,nx+2*haloη,ny+2*haloη)        # tendency of η without time step
            )
    return Tendencies{typeof(NT)}(NT)
end

function VolumeFluxes{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo,haloη = G

    NT = (  h = Array{T,2}(undef,nx+2*haloη,ny+2*haloη),        # layer thickness
            h_u = Array{T,2}(undef,nx+2*haloη-1,ny+2*haloη),    # layer thickness on u-grid
            U = Array{T,2}(undef,nx+2*haloη-1,ny+2*haloη),      # U=uh volume flux

            h_v = Array{T,2}(undef,nx+2*haloη,ny+2*haloη-1),    # layer thickness on v-grid
            V = Array{T,2}(undef,nx+2*haloη,ny+2*haloη-1),      # V=vh volume flux

            dUdx = Array{T,2}(undef,nx+2*haloη-2,ny+2*haloη),   # gradients thereof
            dVdy = Array{T,2}(undef,nx+2*haloη,ny+2*haloη-2)
            )
    return VolumeFluxes{typeof(NT)}(NT)
end

function Vorticity{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo,haloη = G

    NT = (  h_q = Array{T,2}(undef,nx+2*haloη-1,ny+2*haloη-1),  # layer thickness h interpolated on q-grid
            q = Array{T,2}(undef,nx+2*haloη-1,ny+2*haloη-1),    # potential vorticity

            q_v = Array{T,2}(undef,nx+2*haloη-2,ny+2*haloη-1),  # q interpolated on v-grid
            U_v = Array{T,2}(undef,nx+2*haloη-2,ny+2*haloη-1),  # mass flux U=uh on v-grid

            q_u = Array{T,2}(undef,nx+2*haloη-1,ny+2*haloη-2),  # q interpolated on u-grid
            V_u = Array{T,2}(undef,nx+2*haloη-1,ny+2*haloη-2),  # mass flux V=vh on v-grid

            qhu = Array{T,2}(undef,nvx,nvy),            # potential vorticity advection term u-component
            qhv = Array{T,2}(undef,nux,nuy),            # potential vorticity advection term v-component

            u_v = Array{T,2}(undef,nux+2*halo-1,nuy+2*halo-1),  # u-velocity on v-grid
            v_u = Array{T,2}(undef,nvx+2*halo-1,nvy+2*halo-1),  # v-velocity on u-grid

            dudx = Array{T,2}(undef,nux+2*halo-1,nuy+2*halo),   # ∂u/∂x
            dudy = Array{T,2}(undef,nux+2*halo,nuy+2*halo-1),   # ∂u/∂y

            dvdx = Array{T,2}(undef,nvx+2*halo-1,nvy+2*halo),   # ∂v/∂x
            dvdy = Array{T,2}(undef,nvx+2*halo,nvy+2*halo-1)    # ∂v/∂y
            )
    return Vorticity{typeof(NT)}(NT)
end

function Bernoulli{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo,haloη = G

    NT = (  u² = Array{T,2}(undef,nux+2*halo,nuy+2*halo),       # u-velocity squared
            v² = Array{T,2}(undef,nvx+2*halo,nvy+2*halo),       # v-velocity squared

            KEu = Array{T,2}(undef,nux+2*halo-1,nuy+2*halo),    # u-velocity squared on T-grid
            KEv = Array{T,2}(undef,nvx+2*halo,nvy+2*halo-1),    # v-velocity squared on T-grid

            p = Array{T,2}(undef,nx+2*haloη,ny+2*haloη),        # Bernoulli potential
            dpdx = Array{T,2}(undef,nx+2*haloη-1,ny+2*haloη),   # ∂p/∂x
            dpdy = Array{T,2}(undef,nx+2*haloη,ny+2*haloη-1)    # ∂p/∂y
            )
    return Bernoulli{typeof(NT)}(NT)
end

function Bottomdrag{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo, haloη = G

    NT = (  sqrtKE = Array{T,2}(undef,nx+2*haloη,ny+2*haloη),       # sqrt of kinetic energy
            sqrtKE_u = Array{T,2}(undef,nx+2*haloη-1,ny+2*haloη),   # interpolated on u-grid
            sqrtKE_v = Array{T,2}(undef,nx+2*haloη,ny+2*haloη-1),   # interpolated on v-grid

            Bu = Array{T,2}(undef,nx+2*haloη-1,ny+2*haloη),         # bottom friction term u-component
            Bv = Array{T,2}(undef,nx+2*haloη,ny+2*haloη-1)          # bottom friction term v-component
            )
    return Bottomdrag{typeof(NT)}(NT)
end

function ArakawaHsu{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo,haloη = G

    # Linear combination of potential vorticity
    NT = (  qα = Array{T,2}(undef,nx+2*haloη-2,ny+2*haloη-2),
            qβ = Array{T,2}(undef,nx+2*haloη-1,ny+2*haloη-2),
            qγ = Array{T,2}(undef,nx+2*haloη-1,ny+2*haloη-2),
            qδ = Array{T,2}(undef,nx+2*haloη-2,ny+2*haloη-2)
            )
    return ArakawaHsu{typeof(NT)}(NT)
end

function Laplace{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo = G

    NT = (  Lu = Array{T,2}(undef,nux+2*halo-2,nuy+2*halo-2),       # ∇²u
            Lv = Array{T,2}(undef,nvx+2*halo-2,nvy+2*halo-2),       # ∇²v

            # Derivatives of Lu,Lv
            dLudx = Array{T,2}(undef,nux+2*halo-3,nuy+2*halo-2),
            dLudy = Array{T,2}(undef,nux+2*halo-2,nuy+2*halo-3),
            dLvdx = Array{T,2}(undef,nvx+2*halo-3,nvy+2*halo-2),
            dLvdy = Array{T,2}(undef,nvx+2*halo-2,nvy+2*halo-3)
            )
    return Laplace{typeof(NT)}(NT)
end

function Smagorinsky{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo,haloη = G

    NT = (  DT = Array{T,2}(undef,nx+2*haloη,ny+2*haloη),       # Tension squared (on the T-grid)
            DS = Array{T,2}(undef,nx+2*haloη,ny+2*haloη),       # Shearing strain squared (on the T-grid)
            νSmag = Array{T,2}(undef,nx+2*haloη,ny+2*haloη),    # Viscosity coefficient

            # Tension squared on the q-grid
            DS_q = Array{T,2}(undef,nvx+2*halo-1,nvy+2*halo),

            # Smagorinsky viscosity coefficient on the q-grid
            νSmag_q = Array{T,2}(undef,nx+2*haloη-1,ny+2*haloη-1),

            # Entries of the Smagorinsky viscous tensor
            S12 = Array{T,2}(undef,nx+2*haloη-1,ny+2*haloη-1),
            S21 = Array{T,2}(undef,nx+2*haloη-1,ny+2*haloη-1),

            S11 = Array{T,2}(undef,nux+2*halo-3,nuy+2*halo-2),
            S22 = Array{T,2}(undef,nvx+2*halo-2,nvy+2*halo-3),

            # u- and v-components 1 and 2 of the biharmonic diffusion tendencies
            LLu1 = Array{T,2}(undef,nux+2*halo-4,nuy+2*halo-2),
            LLu2 = Array{T,2}(undef,nx+1,ny),

            LLv1 = Array{T,2}(undef,nx,ny+1),
            LLv2 = Array{T,2}(undef,nvx+2*halo-2,nvy+2*halo-4)
            )
    return Smagorinsky{typeof(NT)}(NT)
end

function SemiLagrange{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo,halosstx,halossty = G

    NT = (  xd = Array{T,2}(undef,nx,ny),                       # departure points x-coord
            yd = Array{T,2}(undef,nx,ny),                       # departure points y-coord

            um = Array{T,2}(undef,nux+2*halo,nuy+2*halo),       # u-velocity temporal mid-point
            vm = Array{T,2}(undef,nvx+2*halo,nvy+2*halo),       # v-velocity temporal mid-point

            u_T = Array{T,2}(undef,nux+2*halo-1,nuy+2*halo),    # u-velocity interpolated on T-grid
            um_T = Array{T,2}(undef,nux+2*halo-1,nuy+2*halo),   # um interpolated on T-grid
            v_T = Array{T,2}(undef,nvx+2*halo,nvy+2*halo-1),    # v-velocity interpolated on T-grid
            vm_T = Array{T,2}(undef,nvx+2*halo,nvy+2*halo-1),   # vm interpolated on T-grid

            uinterp = Array{T,2}(undef,nx,ny),                  # u interpolated on mid-point xd,yd
            vinterp = Array{T,2}(undef,nx,ny),                  # v interpolated on mid-point xd,yd

            ssti = Array{T,2}(undef,nx+2*halosstx,ny+2*halossty)# sst interpolated on departure points
            )
    return SemiLagrange{typeof(NT)}(NT)
end

abstract type NamedTupleInStruct end

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

    NT = (  u0 = zeros(T,nux+2*halo,nuy+2*halo),       # u-velocities for RK updates
            u1 = zeros(T,nux+2*halo,nuy+2*halo),

            v0 = zeros(T,nvx+2*halo,nvy+2*halo),       # v-velocities for RK updates
            v1 = zeros(T,nvx+2*halo,nvy+2*halo),

            η0 = zeros(T,nx+2*haloη,ny+2*haloη),       # sea surface height for RK updates
            η1 = zeros(T,nx+2*haloη,ny+2*haloη)
            )
    return RungeKutta{typeof(NT)}(NT)
end

function Tendencies{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo,haloη = G

    NT = (  du = zeros(T,nux+2*halo,nuy+2*halo),       # tendency of u without time step
            dv = zeros(T,nvx+2*halo,nvy+2*halo),       # tendency of v without time step
            dη = zeros(T,nx+2*haloη,ny+2*haloη)        # tendency of η without time step
            )
    return Tendencies{typeof(NT)}(NT)
end

function VolumeFluxes{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo,haloη = G

    NT = (  h = zeros(T,nx+2*haloη,ny+2*haloη),        # layer thickness
            h_u = zeros(T,nx+2*haloη-1,ny+2*haloη),    # layer thickness on u-grid
            U = zeros(T,nx+2*haloη-1,ny+2*haloη),      # U=uh volume flux

            h_v = zeros(T,nx+2*haloη,ny+2*haloη-1),    # layer thickness on v-grid
            V = zeros(T,nx+2*haloη,ny+2*haloη-1),      # V=vh volume flux

            dUdx = zeros(T,nx+2*haloη-2,ny+2*haloη),   # gradients thereof
            dVdy = zeros(T,nx+2*haloη,ny+2*haloη-2)
            )
    return VolumeFluxes{typeof(NT)}(NT)
end

function Vorticity{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo,haloη = G

    NT = (  h_q = zeros(T,nx+2*haloη-1,ny+2*haloη-1),  # layer thickness h interpolated on q-grid
            q = zeros(T,nx+2*haloη-1,ny+2*haloη-1),    # potential vorticity

            q_v = zeros(T,nx+2*haloη-2,ny+2*haloη-1),  # q interpolated on v-grid
            U_v = zeros(T,nx+2*haloη-2,ny+2*haloη-1),  # mass flux U=uh on v-grid

            q_u = zeros(T,nx+2*haloη-1,ny+2*haloη-2),  # q interpolated on u-grid
            V_u = zeros(T,nx+2*haloη-1,ny+2*haloη-2),  # mass flux V=vh on v-grid

            qhu = zeros(T,nvx,nvy),            # potential vorticity advection term u-component
            qhv = zeros(T,nux,nuy),            # potential vorticity advection term v-component

            u_v = zeros(T,nux+2*halo-1,nuy+2*halo-1),  # u-velocity on v-grid
            v_u = zeros(T,nvx+2*halo-1,nvy+2*halo-1),  # v-velocity on u-grid

            dudx = zeros(T,nux+2*halo-1,nuy+2*halo),   # ∂u/∂x
            dudy = zeros(T,nux+2*halo,nuy+2*halo-1),   # ∂u/∂y

            dvdx = zeros(T,nvx+2*halo-1,nvy+2*halo),   # ∂v/∂x
            dvdy = zeros(T,nvx+2*halo,nvy+2*halo-1)    # ∂v/∂y
            )
    return Vorticity{typeof(NT)}(NT)
end

function Bernoulli{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo,haloη = G

    NT = (  u² = zeros(T,nux+2*halo,nuy+2*halo),       # u-velocity squared
            v² = zeros(T,nvx+2*halo,nvy+2*halo),       # v-velocity squared

            KEu = zeros(T,nux+2*halo-1,nuy+2*halo),    # u-velocity squared on T-grid
            KEv = zeros(T,nvx+2*halo,nvy+2*halo-1),    # v-velocity squared on T-grid

            p = zeros(T,nx+2*haloη,ny+2*haloη),        # Bernoulli potential
            dpdx = zeros(T,nx+2*haloη-1,ny+2*haloη),   # ∂p/∂x
            dpdy = zeros(T,nx+2*haloη,ny+2*haloη-1)    # ∂p/∂y
            )
    return Bernoulli{typeof(NT)}(NT)
end

function Bottomdrag{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo, haloη = G

    NT = (  sqrtKE = zeros(T,nx+2*haloη,ny+2*haloη),       # sqrt of kinetic energy
            sqrtKE_u = zeros(T,nx+2*haloη-1,ny+2*haloη),   # interpolated on u-grid
            sqrtKE_v = zeros(T,nx+2*haloη,ny+2*haloη-1),   # interpolated on v-grid

            Bu = zeros(T,nx+2*haloη-1,ny+2*haloη),         # bottom friction term u-component
            Bv = zeros(T,nx+2*haloη,ny+2*haloη-1)          # bottom friction term v-component
            )
    return Bottomdrag{typeof(NT)}(NT)
end

function ArakawaHsu{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo,haloη = G

    # Linear combination of potential vorticity
    NT = (  qα = zeros(T,nx+2*haloη-2,ny+2*haloη-2),
            qβ = zeros(T,nx+2*haloη-1,ny+2*haloη-2),
            qγ = zeros(T,nx+2*haloη-1,ny+2*haloη-2),
            qδ = zeros(T,nx+2*haloη-2,ny+2*haloη-2)
            )
    return ArakawaHsu{typeof(NT)}(NT)
end

function Laplace{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo = G

    NT = (  Lu = zeros(T,nux+2*halo-2,nuy+2*halo-2),       # ∇²u
            Lv = zeros(T,nvx+2*halo-2,nvy+2*halo-2),       # ∇²v

            # Derivatives of Lu,Lv
            dLudx = zeros(T,nux+2*halo-3,nuy+2*halo-2),
            dLudy = zeros(T,nux+2*halo-2,nuy+2*halo-3),
            dLvdx = zeros(T,nvx+2*halo-3,nvy+2*halo-2),
            dLvdy = zeros(T,nvx+2*halo-2,nvy+2*halo-3)
            )
    return Laplace{typeof(NT)}(NT)
end

function Smagorinsky{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo,haloη = G

    NT = (  DT = zeros(T,nx+2*haloη,ny+2*haloη),       # Tension squared (on the T-grid)
            DS = zeros(T,nx+2*haloη,ny+2*haloη),       # Shearing strain squared (on the T-grid)
            νSmag = zeros(T,nx+2*haloη,ny+2*haloη),    # Viscosity coefficient

            # Tension squared on the q-grid
            DS_q = zeros(T,nvx+2*halo-1,nvy+2*halo),

            # Smagorinsky viscosity coefficient on the q-grid
            νSmag_q = zeros(T,nx+2*haloη-1,ny+2*haloη-1),

            # Entries of the Smagorinsky viscous tensor
            S12 = zeros(T,nx+2*haloη-1,ny+2*haloη-1),
            S21 = zeros(T,nx+2*haloη-1,ny+2*haloη-1),

            S11 = zeros(T,nux+2*halo-3,nuy+2*halo-2),
            S22 = zeros(T,nvx+2*halo-2,nvy+2*halo-3),

            # u- and v-components 1 and 2 of the biharmonic diffusion tendencies
            LLu1 = zeros(T,nux+2*halo-4,nuy+2*halo-2),
            LLu2 = zeros(T,nx+1,ny),

            LLv1 = zeros(T,nx,ny+1),
            LLv2 = zeros(T,nvx+2*halo-2,nvy+2*halo-4)
            )
    return Smagorinsky{typeof(NT)}(NT)
end

function SemiLagrange{T}(G::Grid) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = G
    @unpack halo,halosstx,halossty = G

    NT = (  xd = zeros(T,nx,ny),                       # departure points x-coord
            yd = zeros(T,nx,ny),                       # departure points y-coord

            um = zeros(T,nux+2*halo,nuy+2*halo),       # u-velocity temporal mid-point
            vm = zeros(T,nvx+2*halo,nvy+2*halo),       # v-velocity temporal mid-point

            u_T = zeros(T,nux+2*halo-1,nuy+2*halo),    # u-velocity interpolated on T-grid
            um_T = zeros(T,nux+2*halo-1,nuy+2*halo),   # um interpolated on T-grid
            v_T = zeros(T,nvx+2*halo,nvy+2*halo-1),    # v-velocity interpolated on T-grid
            vm_T = zeros(T,nvx+2*halo,nvy+2*halo-1),   # vm interpolated on T-grid

            uinterp = zeros(T,nx,ny),                  # u interpolated on mid-point xd,yd
            vinterp = zeros(T,nx,ny),                  # v interpolated on mid-point xd,yd

            ssti = zeros(T,nx+2*halosstx,ny+2*halossty)# sst interpolated on departure points
            )
    return SemiLagrange{typeof(NT)}(NT)
end

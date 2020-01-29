""" Extends the matrices u,v,η,sst with a halo of ghost points for boundary conditions."""
function add_halo(  u::Array{T,2},
                    v::Array{T,2},
                    η::Array{T,2},
                    sst::Array{T,2},
                    S::ModelSetup) where {T<:AbstractFloat}

    @unpack nx,ny,nux,nuy,nvx,nvy = S.grid
    @unpack halo,haloη,halosstx,halossty = S.grid

    # Add zeros to satisfy kinematic boundary conditions
    u = cat(zeros(T,halo,nuy),u,zeros(T,halo,nuy),dims=1)
    u = cat(zeros(T,nux+2*halo,halo),u,zeros(T,nux+2*halo,halo),dims=2)

    v = cat(zeros(T,halo,nvy),v,zeros(T,halo,nvy),dims=1)
    v = cat(zeros(T,nvx+2*halo,halo),v,zeros(T,nvx+2*halo,halo),dims=2)

    η = cat(zeros(T,haloη,ny),η,zeros(T,haloη,ny),dims=1)
    η = cat(zeros(T,nx+2*haloη,haloη),η,zeros(T,nx+2*haloη,haloη),dims=2)

    sst = cat(zeros(T,halosstx,ny),sst,zeros(T,halosstx,ny),dims=1)
    sst = cat(zeros(T,nx+2*halosstx,halossty),sst,zeros(T,nx+2*halosstx,halossty),dims=2)

    ghost_points!(u,v,η,S)
    ghost_points_sst!(sst,S)
    return u,v,η,sst
end

"""Cut off the halo from the prognostic variables."""
function remove_halo(   u::Array{T,2},
                        v::Array{T,2},
                        η::Array{T,2},
                        sst::Array{T,2},
                        S::ModelSetup) where {T<:AbstractFloat}

    @unpack halo,haloη,halosstx,halossty = S.grid

    uview = u[halo+1:end-halo,halo+1:end-halo]
    vview = v[halo+1:end-halo,halo+1:end-halo]
    ηview = η[haloη+1:end-haloη,haloη+1:end-haloη]
    sstview = sst[halosstx+1:end-halosstx,halossty+1:end-halossty]

    return uview,vview,ηview,sstview
end

""" Copy ghost points for u from inside to the halo in the nonperiodic case. """
function ghost_points_u_nonperiodic!(C::Constants,u::AbstractMatrix)

    @unpack one_minus_α = C

    # tangential boundary condition
    @views @inbounds u[:,1] .= one_minus_α*u[:,4]
    @views @inbounds u[:,2] .= one_minus_α*u[:,3]
    @views @inbounds u[:,end-1] .= one_minus_α*u[:,end-2]
    @views @inbounds u[:,end] .= one_minus_α*u[:,end-3]
end

""" Copy ghost points for u from inside to the halo in the periodic case. """
function ghost_points_u_periodic!(::Type{Tcomm},u::Array{T,2},C::Constants) where {Tcomm,T}

    @unpack one_minus_α = C

    # periodic bc   - mimic MPI communication with (reduced precision) type Tcomm
    @views @inbounds u[1,:] .= T.(Tcomm.(u[end-3,:]))
    @views @inbounds u[2,:] .= T.(Tcomm.(u[end-2,:]))
    @views @inbounds u[end-1,:] .= T.(Tcomm.(u[3,:]))
    @views @inbounds u[end,:] .= T.(Tcomm.(u[4,:]))

    # tangential bc
    @views @inbounds u[:,1] .= one_minus_α*u[:,4]
    @views @inbounds u[:,2] .= one_minus_α*u[:,3]
    @views @inbounds u[:,end-1] .= one_minus_α*u[:,end-2]
    @views @inbounds u[:,end] .= one_minus_α*u[:,end-3]

end

""" Copy ghost points for v from inside to the halo in the nonperiodic case. """
function ghost_points_v_nonperiodic!(C::Constants,v::AbstractMatrix)

    @unpack one_minus_α = C

    # tangential boundary condition
    @views @inbounds v[1,:] .= one_minus_α*v[4,:]
    @views @inbounds v[2,:] .= one_minus_α*v[3,:]
    @views @inbounds v[end-1,:] .= one_minus_α*v[end-2,:]
    @views @inbounds v[end,:] .= one_minus_α*v[end-3,:]
end

""" Copy ghost points for v from inside to the halo in the periodic case. """
function ghost_points_v_periodic!(::Type{Tcomm},v::Array{T,2}) where {Tcomm,T}

    # tangential boundary condition - mimic MPI communication with (reduced precision) type Tcomm
    @views @inbounds v[1,:] .= T.(Tcomm.(v[end-3,:]))
    @views @inbounds v[2,:] .= T.(Tcomm.(v[end-2,:]))
    @views @inbounds v[end-1,:] .= T.(Tcomm.(v[3,:]))
    @views @inbounds v[end,:] .= T.(Tcomm.(v[4,:]))
end

""" Copy ghost points for η from inside to the halo in the nonperiodic case. """
function ghost_points_η_nonperiodic!(η::AbstractMatrix)

    # assume no gradients of η across solid boundaries
    # the 4 corner points are copied twice, but it's faster!
    @views @inbounds η[1,:] .= η[2,:]
    @views @inbounds η[end,:] .= η[end-1,:]

    @views @inbounds η[:,1] .= η[:,2]
    @views @inbounds η[:,end] .= η[:,end-1]
end

""" Copy ghost points for η from inside to the halo in the periodic case. """
function ghost_points_η_periodic!(::Type{Tcomm},η::Array{T,2}) where {Tcomm,T}

    # corner points are copied twice, but it's faster!
    # mimic MPI communication with (reduced precision) type Tcomm
    @views @inbounds η[1,:] .= T.(Tcomm.(η[end-1,:]))
    @views @inbounds η[end,:] .= T.(Tcomm.(η[2,:]))

    @views @inbounds η[:,1] .= T.(Tcomm.(η[:,2]))
    @views @inbounds η[:,end] .= T.(Tcomm.(η[:,end-1]))
end

""" Copy ghost points for η from inside to the halo in the nonperiodic case. """
function ghost_points_sst_nonperiodic!(G::Grid,sst::AbstractMatrix)

    @unpack halosstx,halossty = G

    # assume no gradients of η across solid boundaries
    # the 4 corner points are copied twice, but it's faster!
    for i ∈ 1:halosstx
        @views @inbounds sst[i,:] .= sst[halosstx+1,:]
        @views @inbounds sst[end-i+1,:] .= sst[end-halosstx,:]
    end

    for j ∈ 1:halossty
        @views @inbounds sst[:,j] .= sst[:,halossty+1]
        @views @inbounds sst[:,end-j+1] .= sst[:,end-halossty]
    end
end

""" Copy ghost points for η from inside to the halo in the periodic case. """
function ghost_points_sst_periodic!(::Type{Tcomm},sst::Array{T,2},G::Grid) where {Tcomm,T}

    @unpack halosstx,halossty = G

    # corner points are copied twice, but it's faster!
    for i ∈ 1:halosstx
        @views @inbounds sst[i,:] .= T.(Tcomm.(sst[end-2*halosstx+i,:]))
        @views @inbounds sst[end-halosstx+i,:] .= T.(Tcomm.(sst[halosstx+i,:]))
    end

    for j ∈ 1:halossty
        @views @inbounds sst[:,j] .= sst[:,halossty+1]
        @views @inbounds sst[:,end-j+1] .= sst[:,end-halossty]
    end
end

"""Decide on boundary condition P.bc which ghost point function to execute."""
function ghost_points!( u::AbstractMatrix,
                        v::AbstractMatrix,
                        η::AbstractMatrix,
                        S::ModelSetup)

    @unpack bc = S.parameters
    C = S.constants

    if bc == "periodic"
        @unpack Tcomm = S.parameters
        ghost_points_u_periodic!(Tcomm,u,C)
        ghost_points_v_periodic!(Tcomm,v)
        ghost_points_η_periodic!(Tcomm,η)
    else
        ghost_points_u_nonperiodic!(C,u)
        ghost_points_v_nonperiodic!(C,v)
        ghost_points_η_nonperiodic!(η)
    end
end

"""Decide on boundary condition P.bc which ghost point function to execute."""
function ghost_points_sst!(sst::AbstractMatrix,S::ModelSetup)

    @unpack bc = S.parameters
    G = S.grid

    if bc == "periodic"
        @unpack Tcomm = S.parameters
        ghost_points_sst_periodic!(Tcomm,sst,G)
    else
        ghost_points_sst_nonperiodic!(G,sst)
    end
end

"""
This function will compute an additional forcing term S that will act
as an eddy parameterization. This is done by
(1) Adding a new structure that will preallocate all of the operators I need (called ZB_momentum)
(2) Computing the main operators that appear in the new term: ζ, D, Dhat. ζ is the relative vorticity,
D is the shear deformation of the flow field (both of these live on cell corners), and Dhat is the stretch
deformation of the flow field (this lives on the cell centers). The computation of these operators was 
extensively checked for accuracy by MK and SW
(3) Using ζ, D, and Dhat, then compute the individual pieces of S that appear, comments
below explain on which grid they live and so on

There is also an option to apply a filter to S. The first few inline functions relate to this 
operation. The application of the filter (or convolutional kernal) is done in three stages: 
(1) applying the filter only to the internal points
(2) applying the filter to points on the boundary, not including corners
(3) applying the filter to points on the corners

There are two options about where to add this function: it can be added during each RK step with rhs_nonlinear!, or it can be added once during the 
diffusion computation with add_drag_diff_tendencies!. The option is decided when default_parameters is set.
"""

@inline function Gconvolve(array)
    return (array[1, 1] + 2*array[1, 2] + array[1, 3] + 2*array[2, 1] + 4*array[2, 2] + 2*array[2, 3] + array[3, 1] + 2*array[3, 2] + array[3, 3]) / 16
end
@inline function Gconvolve12all(array)
    return (array[1, 1] + 2*array[1, 2] + array[1, 3] + 2*array[2, 1] + 4*array[2, 2] + 2*array[2, 3]) / 12
end
@inline function Gconvolve23all(array)
    return (2*array[1, 1] + 4*array[1, 2] + 2*array[1, 3] + array[2, 1] + 2*array[2, 2] + array[2, 3]) / 12
end
@inline function Gconvolveall12(array)
    return (array[1, 1] + 2*array[1, 2] + 2*array[2, 1] + 4*array[2, 2] + array[3, 1] + 2*array[3, 2]) / 12
end
@inline function Gconvolveall23(array)
    return (2*array[1, 1] + array[1, 2] + 4*array[2, 1] + 2*array[2, 2] + 2*array[3, 1] + array[3, 2]) / 12
end

function apply_filter(N,ζsqT,ζDT,ζDhat,ζD_filtered,ζDhat_filtered,trace_filtered)

    mT,nT = size(ζsqT)
    mq,nq = size(ζDhat)

    for i = 1:N

        @inbounds for j ∈ 2:nT-1 
            for k ∈ 2:mT-1 
                ζsqT[k,j] = Gconvolve(@view ζsqT[k-1:k+1,j-1:j+1])
                ζDT[k,j] = Gconvolve(@view ζDT[k-1:k+1,j-1:j+1])
            end
        end

        @inbounds for h ∈ 2:nT-1
            ζsqT[1,h] = Gconvolve23all(@view ζsqT[1:2,h-1:h+1])
            ζsqT[mT,h] = Gconvolve12all(@view ζsqT[mT-1:mT,h-1:h+1])
            ζDT[1,h] = Gconvolve23all(@view ζDT[1:2,h-1:h+1])
            ζDT[mT,h] = Gconvolve12all(@view ζDT[mT-1:mT,h-1:h+1])
        end

        @inbounds for v ∈ 2:mT-1
            ζsqT[v,1] = Gconvolveall23(@view ζsqT[v-1:v+1,1:2])
            ζsqT[v,nT] = Gconvolveall12(@view ζsqT[v-1:v+1,nT-1:nT])
            ζDT[v,1] = Gconvolveall23(@view ζDT[v-1:v+1,1:2])
            ζDT[v,nT] = Gconvolveall12(@view ζDT[v-1:v+1,nT-1:nT])
        end

        ζsqT[1,1] = (4*ζsqT[1,1] + 2*ζsqT[1,2] + 2*ζsqT[2,1] + ζsqT[2,2])/8
        ζsqT[1,nT] = (2*ζsqT[1,nT-1] + 4*ζsqT[1,nT] + ζsqT[2,nT-1] + 2*ζsqT[2,nT])/8
        ζsqT[mT,1] = (2*ζsqT[mT-1,1] + ζsqT[mT-1,2] + 4*ζsqT[mT,1] + 2*ζsqT[mT,2])/8
        ζsqT[mT,nT] = (ζsqT[mT-1,nT-1] + 2*ζsqT[mT-1,nT] + 2*ζsqT[mT,nT-1] + 4*ζsqT[mT,nT])/8

        ζDT[1,1] = (4*ζDT[1,1] + 2*ζDT[1,2] + 2*ζDT[2,1] + ζDT[2,2])/8
        ζDT[1,nT] = (2*ζDT[1,nT-1] + 4*ζDT[1,nT] + ζDT[2,nT-1] + 2*ζDT[2,nT])/8
        ζDT[mT,1] = (2*ζDT[mT-1,1] + ζDT[mT-1,2] + 4*ζDT[mT,1] + 2*ζDT[mT,2])/8
        ζDT[mT,nT] = (ζDT[mT-1,nT-1] + 2*ζDT[mT-1,nT] + 2*ζDT[mT,nT-1] + 4*ζDT[mT,nT])/8

        @inbounds for j ∈ 2:nq-1
            for k ∈ 2:mq-1
                ζDhat[k,j] = Gconvolve(@view ζDhat[k-1:k+1,j-1:j+1])
            end
        end

        @inbounds for h ∈ 2:nq-1
            ζDhat[1,h] = Gconvolve23all(@view ζDhat[1:2,h-1:h+1])
            ζDhat[mq,h] = Gconvolve12all(@view ζDhat[mq-1:mq,h-1:h+1])
        end

        @inbounds for v ∈ 2:mq-1
            ζDhat[v,1] = Gconvolveall23(@view ζDhat[v-1:v+1,1:2])
            ζDhat[v,nq] = Gconvolveall12(@view ζDhat[v-1:v+1,nq-1:nq])
        end

        ζDhat[1,1] = (4*ζDhat[1,1] + 2*ζDhat[1,2] + 2*ζDhat[2,1] + ζDhat[2,2])/8
        ζDhat[1,nq] = (2*ζDhat[1,nq-1] + 4*ζDhat[1,nq] + ζDhat[2,nq-1] + 2*ζDhat[2,nq])/8
        ζDhat[mq,1] = (2*ζDhat[mq-1,1] + ζDhat[mq-1,2] + 4*ζDhat[mq,1] + 2*ζDhat[mq,2])/8
        ζDhat[mq,nq] = (ζDhat[mq-1,nq-1] + 2*ζDhat[mq-1,nq] + 2*ζDhat[mq,nq-1] + 4*ζDhat[mq,nq])/8

    end

    trace_filtered .= ζsqT
    ζD_filtered .= ζDT
    ζDhat_filtered .= ζDhat

    return nothing

end

function ZB_forcing(u, v, S, Diag)

    @unpack zb_filtered, N  = S.parameters
    @unpack γ₀, ζ, ζsq, D, Dsq, Dhat, Dhatsq, Dhatq = Diag.ZBVars
    @unpack ζD, ζDT, ζDhat, ζsqT, trace = Diag.ZBVars
    @unpack ζpDT = Diag.ZBVars
    @unpack dudx, dudy, dvdx, dvdy = Diag.ZBVars

    @unpack dζDdx, dζDhatdy, dtracedx = Diag.ZBVars
    @unpack dζDhatdx, dζDdy, dtracedy = Diag.ZBVars
    @unpack S_u, S_v = Diag.ZBVars
    @unpack Δ, scale, f₀ = S.grid

    @unpack Ker = Diag.ZBVars
    @unpack ζD_filtered, ζDhat_filtered, trace_filtered = Diag.ZBVars

    @unpack halo, haloη, ep, nux, nuy, nvx, nvy = S.grid

    κ_BC = - γ₀ * Δ^2

    ∂x!(dudx, u)
    ∂y!(dudy, u)

    ∂x!(dvdx, v)
    ∂y!(dvdy, v)

    mq,nq = size(ζ)
    mTh,nTh = size(Dhat)
    mT,nT = size(trace)

    ##### CHECK #########
    # @boundscheck (mq+2,nq+2) == size(dvdx) || throw(BoundsError())
    # @boundscheck (mq+2+ep,nq+2) == size(dudy) || throw(BoundsError())

    # Relative vorticity and shear deformation, cell corners
    @inbounds for j ∈ 1:nq
        for k ∈ 1:mq
            ζ[k,j] = dvdx[k+1,j+1] - dudy[k+1,j+1]
            D[k,j] = dudy[k+1,j+1] + dvdx[k+1,j+1]
        end
    end

    # Stretch deformation, cell centers (with halo)
    @inbounds for j ∈ 1:nTh
        for k ∈ 1:mTh
            Dhat[k,j] = dudx[k,j+1] - dvdy[k+1,j]
        end
    end

    ζsq .= κ_BC .* ζ.^2
    Dsq .= D.^2
    Dhatsq .= Dhat.^2

    # Trace computation (second term in forcing term), only keeping ζ^2 and no other terms
    Ixy!(ζsqT, ζsq)

    # Computing ζ ⋅ D and placing on cell centers
    ζD .= κ_BC .* (ζ .* D)
    Ixy!(ζDT, ζD)

    # Computing ζ ⋅ Dhat, cell corners
    Ixy!(Dhatq, Dhat)
    for kj in eachindex(ζDhat,ζ,Dhatq)
        @inbounds ζDhat[kj] = κ_BC * ζ[kj] * Dhatq[kj]
    end

    if zb_filtered

        apply_filter(N, ζsqT, ζDT, ζDhat, ζD_filtered, ζDhat_filtered, trace_filtered)

        ∂x!(dζDdx, ζD_filtered)
        ∂y!(dζDhatdy, ζDhat_filtered)
        ∂x!(dtracedx, trace_filtered)

        ∂x!(dζDhatdx, ζDhat_filtered)
        ∂y!(dζDdy, ζD_filtered)
        ∂y!(dtracedy, trace_filtered)

    else

        ∂x!(dζDdx, ζDT)
        ∂y!(dζDhatdy, ζDhat)
        ∂x!(dtracedx, ζsqT)

        ∂x!(dζDhatdx, ζDhat)
        ∂y!(dζDdy, ζDT)
        ∂y!(dtracedy, ζsqT)

    end

    s = Δ^2 * scale
    @inbounds for j ∈ 1:nuy
        for k ∈ 1:nux
            S_u[k,j] = (-dζDdx[k,j] + dζDhatdy[k+1,j] + dtracedx[k,j]) / s
        end
    end

    @inbounds for j ∈ 1:nvy
        for k ∈ 1:nvx
            S_v[k,j] = (dζDhatdx[k,j+1] + dζDdy[k,j] + dtracedy[k,j]) / s
        end
    end

end
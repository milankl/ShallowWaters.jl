"""Integrates Juls forward in time."""
function time_integration(  ::Type{T},
                            Prog::PrognosticVars,
                            Diag::DiagnosticVars,
                            S::ModelSetup) where {T<:AbstractFloat}

    @unpack u,v,η,sst = Prog
    @unpack u0,v0,η0 = Diag.RungeKutta
    @unpack u1,v1,η1 = Diag.RungeKutta
    @unpack du,dv,dη = Diag.Tendencies
    @unpack um,vm = Diag.SemiLagrange

    @unpack dynamics,RKo,tracer_advection = S.parameters
    @unpack RKaΔt,RKbΔt = S.constants

    @unpack nt,dtint = S.grid
    @unpack nstep_advcor,nstep_diff,nadvstep,nadvstep_half = S.grid

    if dynamics == "linear"
        Ix!(Diag.VolumeFluxes.h_u,S.forcing.H)
        Iy!(Diag.VolumeFluxes.h_v,S.forcing.H)
    end

    # propagate initial conditions
    copyto!(u0,u)
    copyto!(v0,v)
    copyto!(η0,η)

    # feedback and output
    feedback = feedback_init(S)
    netCDFfiles = NcFiles(Prog,Diag,S)

    nans_detected = false
    t = 0           # model time
    for i = 1:nt

        # ghost point copy for boundary conditions
        ghost_points!(u,v,η,S)
        copyto!(u1,u)
        copyto!(v1,v)
        copyto!(η1,η)

        # Runge-Kutta 4th order / 3rd order
        for rki = 1:RKo
            if rki > 1
                ghost_points!(u1,v1,η1,S)
            end

            rhs!(u1,v1,η1,Diag,S)

            if rki < RKo
                caxb!(u1,u,RKbΔt[rki],du)   #u1 .= u .+ RKb[rki]*Δt*du
                caxb!(v1,v,RKbΔt[rki],dv)   #v1 .= v .+ RKb[rki]*Δt*dv
                caxb!(η1,η,RKbΔt[rki],dη)   #η1 .= η .+ RKb[rki]*Δt*dη
            end

            # sum RK-substeps on the go
            axb!(u0,RKaΔt[rki],du)          #u0 .+= RKa[rki]*Δt*du
            axb!(v0,RKaΔt[rki],dv)          #v0 .+= RKa[rki]*Δt*dv
            axb!(η0,RKaΔt[rki],dη)          #η0 .+= RKa[rki]*Δt*dη
        end

        ghost_points!(u0,v0,η0,S)

        # ADVECTION and CORIOLIS TERMS
        # although included in the tendency of every RK substep,
        # only update every nstep_advcor steps if nstep_advcor > 0
        if dynamics == "nonlinear" &&
            nstep_advcor > 0 &&
            (i % nstep_advcor) == 0
            advection_coriolis!(u0,v0,η0,Diag,S)
        end

        # DIFFUSIVE TERMS - SEMI-IMPLICIT EULER
        # use u0 = u^(n+1) to evaluate tendencies, add to u0 = u^n + rhs
        # evaluate only every nstep_diff time steps
        if (i % nstep_diff) == 0
            bottom_drag!(u0,v0,η0,Diag,S)
            diffusion!(u0,v0,Diag,S)
            add_drag_diff_tendencies!(u0,v0,Diag,S)
        end

        # RK3/4 copy back from substeps
        copyto!(u,u0)
        copyto!(v,v0)
        copyto!(η,η0)
        t += dtint

        # TRACER ADVECTION
        tracer!(i,Prog,Diag,S)

        # feedback and output
        feedback.i = i
        feedback!(Prog,feedback)
        output_nc!(i,netCDFfiles,Prog,Diag,S)

        if feedback.nans_detected
            break
        end
    end

    # finalise feedback and output
    feedback_end!(feedback)
    output_close(netCDFfiles,feedback,S)

    return PrognosticVars{T}(remove_halo(u,v,η,sst,S)...)
end

"""Add to a x multiplied with b. a += x*b """
function axb!(a::Array{T,2},x::T,b::Array{T,2}) where {T<:AbstractFloat}
    m,n = size(a)
    @boundscheck (m,n) == size(b) || throw(BoundsError())

    #TODO @simd?
    @inbounds for j ∈ 1:n
        for i ∈ 1:m
           a[i,j] += x*b[i,j]
        end
    end
end

"""c equals add a to x multiplied with b. c = a + x*b """
function caxb!(c::Array{T,2},a::Array{T,2},x::T,b::Array{T,2}) where {T<:AbstractFloat}
    m,n = size(a)
    @boundscheck (m,n) == size(b) || throw(BoundsError())
    @boundscheck (m,n) == size(c) || throw(BoundsError())

    #TODO @simd?
    @inbounds for j ∈ 1:n
        for i ∈ 1:m
           c[i,j] = a[i,j] + x*b[i,j]
        end
    end
end

"""Integrates ShallowWaters forward in time."""
function time_integration(  Prog::PrognosticVars{Tprog},
                            Diag::DiagnosticVars{T,Tprog},
                            S::ModelSetup{T,Tprog}) where {T<:AbstractFloat,Tprog<:AbstractFloat}

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

    # feedback, output initialisation and storing initial conditions
    feedback = feedback_init(S)
    netCDFfiles = NcFiles(feedback,S)
    output_nc!(0,netCDFfiles,Prog,Diag,S)

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

            # type conversion for mixed precision
            u1rhs = convert(Diag.PrognosticVarsRHS.u,u1)
            v1rhs = convert(Diag.PrognosticVarsRHS.v,v1)
            η1rhs = convert(Diag.PrognosticVarsRHS.η,η1)

            rhs!(u1rhs,v1rhs,η1rhs,Diag,S,t)

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

        # type conversion for mixed precision
        u0rhs = convert(Diag.PrognosticVarsRHS.u,u0)
        v0rhs = convert(Diag.PrognosticVarsRHS.v,v0)
        η0rhs = convert(Diag.PrognosticVarsRHS.η,η0)

        # ADVECTION and CORIOLIS TERMS
        # although included in the tendency of every RK substep,
        # only update every nstep_advcor steps if nstep_advcor > 0
        if dynamics == "nonlinear" && nstep_advcor > 0 && (i % nstep_advcor) == 0
            advection_coriolis!(u0rhs,v0rhs,η0rhs,Diag,S)
        end

        # DIFFUSIVE TERMS - SEMI-IMPLICIT EULER
        # use u0 = u^(n+1) to evaluate tendencies, add to u0 = u^n + rhs
        # evaluate only every nstep_diff time steps
        if (i % nstep_diff) == 0
            bottom_drag!(u0rhs,v0rhs,η0rhs,Diag,S)
            diffusion!(u0rhs,v0rhs,Diag,S)
            add_drag_diff_tendencies!(u0,v0,Diag,S)
        end

        # RK3/4 copy back from substeps
        copyto!(u,u0)
        copyto!(v,v0)
        copyto!(η,η0)
        t += dtint

        # TRACER ADVECTION
        tracer!(i,u0rhs,v0rhs,Prog,Diag,S)

        # feedback and output
        feedback.i = i
        feedback!(Prog,feedback,S)
        output_nc!(i,netCDFfiles,Prog,Diag,S)

        if feedback.nans_detected
            break
        end
    end

    # finalise feedback and output
    feedback_end!(feedback)
    output_close!(netCDFfiles,feedback,S)

    return PrognosticVars{Tprog}(remove_halo(u,v,η,sst,S)...)
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

"""Convert function for two arrays, X1, X2, in case their eltypes differ.
Convert every element from X1 and store it in X2."""
function Base.convert(X2::Array{T2,N},X1::Array{T1,N}) where {T1,T2,N}

    @boundscheck size(X2) == size(X1) || throw(BoundsError())

    @inbounds for i in eachindex(X1)
            X2[i] = T2(X1[i])
    end

    return X2
end


"""Convert function for two arrays, X1, X2, in case their eltypes are identical.
Just pass X1, such that X2 is pointed to the same place in memory."""
function Base.convert(X2::Array{T,N},X1::Array{T,N}) where {T,N}
    @boundscheck size(X2) == size(X1) || throw(BoundsError())
    return X1
end

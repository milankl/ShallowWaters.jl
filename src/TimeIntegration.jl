"""Integrates ShallowWaters forward in time."""
function time_integration(  Prog::PrognosticVars{Tprog},
                            Diag::DiagnosticVars{T,Tprog},
                            S::ModelSetup{T,Tprog}) where {T<:AbstractFloat,Tprog<:AbstractFloat}

    @unpack u,v,η,sst = Prog
    @unpack u0,v0,η0 = Diag.RungeKutta
    @unpack u1,v1,η1 = Diag.RungeKutta
    @unpack du,dv,dη = Diag.Tendencies
    @unpack um,vm = Diag.SemiLagrange

    @unpack dynamics,RKo,RKs,tracer_advection = S.parameters
    @unpack time_scheme = S.parameters
    @unpack RKaΔt,RKbΔt = S.constants
    @unpack Δt_Δ,Δt_Δs = S.constants

    @unpack nt,dtint = S.grid
    @unpack nstep_advcor,nstep_diff,nadvstep,nadvstep_half = S.grid

    # calculate layer thicknesses for initial conditions
    thickness!(Diag.VolumeFluxes.h,η,S.forcing.H)
    Ix!(Diag.VolumeFluxes.h_u,Diag.VolumeFluxes.h)
    Iy!(Diag.VolumeFluxes.h_v,Diag.VolumeFluxes.h)
    Ixy!(Diag.Vorticity.h_q,Diag.VolumeFluxes.h)

    # calculate PV terms for initial conditions
    urhs = convert(Diag.PrognosticVarsRHS.u,u)
    vrhs = convert(Diag.PrognosticVarsRHS.v,v)
    ηrhs = convert(Diag.PrognosticVarsRHS.η,η)
    advection_coriolis!(urhs,vrhs,ηrhs,Diag,S)
    PVadvection!(Diag,S)

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

        if time_scheme == "RK"   # classic RK4,3 or 2

            for rki = 1:RKo
                if rki > 1
                    ghost_points!(u1,v1,η1,S)
                end

                # type conversion for mixed precision
                u1rhs = convert(Diag.PrognosticVarsRHS.u,u1)
                v1rhs = convert(Diag.PrognosticVarsRHS.v,v1)
                η1rhs = convert(Diag.PrognosticVarsRHS.η,η1)

                rhs!(u1rhs,v1rhs,η1rhs,Diag,S,t)          # momentum only
                continuity!(u1rhs,v1rhs,η1rhs,Diag,S,t)   # continuity equation 

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

        elseif time_scheme == "SSPRK2"  # s-stage 2nd order SSPRK

            for rki = 1:RKs
                if rki > 1
                    # TODO technically the ghost point copy for u1,v1 is redundant as done further down
                    ghost_points!(u1,v1,η1,S)
                end

                # type conversion for mixed precision
                u1rhs = convert(Diag.PrognosticVarsRHS.u,u1)
                v1rhs = convert(Diag.PrognosticVarsRHS.v,v1)
                η1rhs = convert(Diag.PrognosticVarsRHS.η,η1)

                rhs!(u1rhs,v1rhs,η1rhs,Diag,S,t)

                # the update step
                axb!(u1,Δt_Δs,du)       # u1 = u1 + Δt/(s-1)*RHS(u1)
                axb!(v1,Δt_Δs,dv)

                # semi-implicit for continuity equation, use u1,v1 to calcualte dη
                ghost_points!(u1,v1,S)
                u1rhs = convert(Diag.PrognosticVarsRHS.u,u1)
                v1rhs = convert(Diag.PrognosticVarsRHS.v,v1)
                continuity!(u1rhs,v1rhs,η1rhs,Diag,S,t)
                axb!(η1,Δt_Δs,dη)       # η1 = η1 + Δt/(s-1)*RHS(u1)
            end

            a = 1/RKs
            b = (RKs-1)/RKs
            cxayb!(u0,a,u,b,u1)
            cxayb!(v0,a,v,b,v1)
            cxayb!(η0,a,η,b,η1)
            
        elseif time_scheme == "SSPRK3"   # 4-stage SSPRK3
        
            for rki = 1:4
                if rki > 1
                    ghost_points!(u1,v1,η1,S)
                end

                # type conversion for mixed precision
                u1rhs = convert(Diag.PrognosticVarsRHS.u,u1)
                v1rhs = convert(Diag.PrognosticVarsRHS.v,v1)
                η1rhs = convert(Diag.PrognosticVarsRHS.η,η1)

                rhs!(u1rhs,v1rhs,η1rhs,Diag,S,t)

                caxb!(u0,u1,Δt_Δ,du)        # store Euler update into u0,v0
                caxb!(v0,v1,Δt_Δ,dv)
                cxab!(u1,1/2,u1,u0)         # average u0,u1 and store in u1
                cxab!(v1,1/2,v1,v0)         # same

                # semi-implicit for continuity equation, use u1,v1 to calcualte dη
                ghost_points!(u1,v1,S)
                u1rhs = convert(Diag.PrognosticVarsRHS.u,u1)
                v1rhs = convert(Diag.PrognosticVarsRHS.v,v1)
                continuity!(u1rhs,v1rhs,η1rhs,Diag,S,t)
                
                caxb!(η0,η1,Δt_Δ,dη)    # store Euler update into η0
                cxab!(η1,1/2,η1,η0)         # average η0,η1 and store in η1

                if rki == 3
                    cxayb!(u1,2/3,u,1/3,u1)
                    cxayb!(v1,2/3,v,1/3,v1)
                    cxayb!(η1,2/3,η,1/3,η1)
                elseif rki == 4
                    copyto!(u0,u1)
                    copyto!(v0,v1)
                    copyto!(η0,η1)
                end
            end
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
            UVfluxes!(u0rhs,v0rhs,η0rhs,Diag,S)
            advection_coriolis!(u0rhs,v0rhs,η0rhs,Diag,S)
        end

        # DIFFUSIVE TERMS - SEMI-IMPLICIT EULER
        # use u0 = u^(n+1) to evaluate tendencies, add to u0 = u^n + rhs
        # evaluate only every nstep_diff time steps
        if (i % nstep_diff) == 0
            bottom_drag!(u0rhs,v0rhs,η0rhs,Diag,S)
            diffusion!(u0rhs,v0rhs,Diag,S)
            add_drag_diff_tendencies!(u0,v0,Diag,S)
            ghost_points!(u0,v0,S)
        end

        t += dtint

        # TRACER ADVECTION
        u0rhs = convert(Diag.PrognosticVarsRHS.u,u0)  # copy back as add_drag_diff_tendencies changed u0,v0
        v0rhs = convert(Diag.PrognosticVarsRHS.v,v0)
        tracer!(i,u0rhs,v0rhs,Prog,Diag,S)

        # feedback and output
        feedback.i = i
        feedback!(Prog,feedback,S)
        output_nc!(i,netCDFfiles,Prog,Diag,S)       # uses u0,v0,η0

        if feedback.nans_detected
            break
        end

        # Copy back from substeps
        copyto!(u,u0)
        copyto!(v,v0)
        copyto!(η,η0)
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

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
           c[i,j] = a[i,j] + x*b[i,j]
        end
    end
end

"""c equals add x multiplied to a plus b. c = x*(a+b) """
function cxab!(c::Array{T,2},x::Real,a::Array{T,2},b::Array{T,2}) where {T<:AbstractFloat}
    m,n = size(a)
    @boundscheck (m,n) == size(b) || throw(BoundsError())
    @boundscheck (m,n) == size(c) || throw(BoundsError())

    xT = T(x)
    @inbounds for j ∈ 1:n
        for i ∈ 1:m
           c[i,j] = xT*(a[i,j] + b[i,j])
        end
    end
end

"""c equals add x multiplied to a plus b. c = x*(a+b) """
function cxayb!(c::Array{T,2},x::Real,a::Array{T,2},y::Real,b::Array{T,2}) where {T<:AbstractFloat}
    m,n = size(a)
    @boundscheck (m,n) == size(b) || throw(BoundsError())
    @boundscheck (m,n) == size(c) || throw(BoundsError())

    xT = T(x)
    yT = T(y)
    @inbounds for j ∈ 1:n
        for i ∈ 1:m
           c[i,j] = xT*a[i,j] + yT*b[i,j]
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

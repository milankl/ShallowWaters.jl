"""Integrates Juls forward in time."""
function time_integration!( Prog::PrognosticVars,
                            Diag::DiagnosticVars,
                            S::ModelSetup)

    @unpack u,v,η,sst = Prog
    @unpack u0,v0,η0 = Diag.RungeKutta
    @unpack u1,v1,η1 = Diag.RungeKutta
    @unpack du,dv,dη = Diag.Tendencies
    @unpack um,vm = Diag.SemiLagrange

    @unpack dynamics,RKo,tracer_advection = S.parameter
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
    #ncs_progn,ncs_tend,ncs_diagn,iout = output_ini(u,v,η,sst,du,dv,dη,qhv,qhu,dpdx,dpdy,dUdx,dVdy,Bu,Bv,LLu1,LLu2,LLv1,LLv2,
    #                                                q,p,dudx,dvdy,dudy,dvdx,Lu,Lv,xd,yd,f_q)

    nans_detected = false
    t = 0           # model time
    for i = 1:nt

        # ghost point copy for boundary conditions
        ghost_points!(P,C,u,v,η)
        copyto!(u1,u)
        copyto!(v1,v)
        copyto!(η1,η)

        # Runge-Kutta 4th order / 3rd order
        for rki = 1:RKo
            if rki > 1
                ghost_points!(P,C,u1,v1,η1)
            end

            rhs!(u1,v1,η1,P,C,G,Diag,Forc)

            if rki < RKo
                caxb!(u1,u,RKbΔt[rki],du)  #u1 .= u .+ RKb[rki]*Δt*du
                caxb!(v1,v,RKbΔt[rki],dv)  #v1 .= v .+ RKb[rki]*Δt*dv
                caxb!(η1,η,RKbΔt[rki],dη)  #η1 .= η .+ RKb[rki]*Δt*dη
            end

            # sum RK-substeps on the go
            axb!(u0,RKaΔt[rki],du)  #u0 .+= RKa[rki]*Δt*du
            axb!(v0,RKaΔt[rki],dv)  #v0 .+= RKa[rki]*Δt*dv
            axb!(η0,RKaΔt[rki],dη)  #η0 .+= RKa[rki]*Δt*dη
        end

        ghost_points!(P,C,u0,v0,η0)

        # ADVECTION and CORIOLIS TERMS
        # although included in the tendency of every RK substep,
        # only update every nstep_advcor steps if nstep_advcor > 0
        if dynamics == "nonlinear" && nstep_advcor > 0 && (i % nstep_advcor) == 0
            advection_coriolis!(u0,v0,η0,P,G,Diag,Forc)
        end

        # DIFFUSIVE TERMS - SEMI-IMPLICIT EULER
        # use u0 = u^(n+1) to evaluate tendencies, add to u0 = u^n + rhs
        # evaluate only every nstep_diff time steps
        if (i % nstep_diff) == 0
            bottom_drag!(u0,v0,η0,P,C,G,Diag,Forc)
            diffusion!(u0,v0,P,C,G,Diag)
            add_drag_diff_tendencies!(u0,v0,G,Diag)
        end

        # RK3/4 copy back from substeps
        copyto!(u,u0)
        copyto!(v,v0)
        copyto!(η,η0)
        t += dtint


        # TRACER ADVECTION
        # mid point (in time) velocity for the advective time step
        # if tracer_advection && ((i+nadvstep_half) % nadvstep) == 0
        #     um .= u
        #     vm .= v
        # end
        #
        # if tracer_advection && (i % nadvstep) == 0
        #     departure!(u,v,u_T,v_T,um,vm,um_T,vm_T,uinterp,vinterp,xd,yd)
        #     adv_sst!(ssti,sst,xd,yd)
        #     if tracer_relaxation
        #         tracer_relax!(ssti,sst_ref,SSTγ)
        #     end
        #     if tracer_consumption
        #         tracer_consumption!(ssti)
        #     end
        #     ghost_points_sst!(ssti)
        #     sst .= ssti
        #
        #     # conserved?
        #     #println(mean(sst[halosstx+1:end-halosstx,halossty+1:end-halossty].*h[haloη+1:end-haloη,haloη+1:end-haloη]))
        # end

        # feedback and output
        t0,nans_detected = feedback(u,v,η,sst,i,t0,nt,nans_detected,progrtxt,P)
        feedback!(Prog,feedback)

        #ncs_diagn = output_diagn_nc(ncs_diagn,i,iout,q,p,dudx,dvdy,dudy,dvdx,Lu,Lv,xd,yd,f_q)
        #ncs_tend = output_tend_nc(ncs_tend,i,iout,du,dv,dη,qhv,qhu,dpdx,dpdy,dUdx,dVdy,Bu,Bv,LLu1,LLu2,LLv1,LLv2)
        #ncs_progn,iout = output_progn_nc(ncs_progn,i,iout,u,v,η,sst)

        if feedback.nans_detected
            break
        end
    end

    # finalise feedback and output
    feedback_end!(feedback)
    #output_close(ncs_progn,ncs_tend,ncs_diagn,progrtxt)
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

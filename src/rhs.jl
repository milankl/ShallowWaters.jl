function rhs!(du,dv,dη,u,v,η,Fx,f_q,H,
             dudx,dvdy,dvdx,dudy,dpdx,dpdy,
             p,KEu,KEv,dUdx,dVdy,
             h,h_u,h_v,h_q,U,V,U_v,V_u,
             adv_u,adv_v,q,q_u,q_v,
             Lu1,Lu2,Lv1,Lv2,dLudx,dLudy,dLvdx,dLvdy,
             shear,νSmag,νSmag_q)

    @views h .= η .+ H
    Ix!(h_u,h)
    Iy!(h_v,h)
    Ixy!(h_q,h)

    Uflux!(U,u,h_u)
    Vflux!(V,v,h_v)

    ∂x!(dudx,u)
    ∂y!(dvdy,v)
    ∂x!(dvdx,v)
    ∂y!(dudy,u)

    ∂x!(dUdx,U)
    ∂y!(dVdy,V)

    # Bernoulli potential
    @views u² .= u.^2
    @views v² .= v.^2
    Ix!(KEu,u²)
    Iy!(KEv,v²)
    Bernoulli!(p,KEu,KEv,η)
    ∂x!(dpdx,p)
    ∂y!(dpdy,p)

    # Potential vorticity
    PV!(q,f_q,dvdx,dudy,h_q)

    # Sadourny, 1975 enstrophy conserving scheme
    Iy!(q_u,q)
    Ix!(q_v,q)
    Ixy!(V_u,V)
    Ixy!(U_v,U)
    PV_adv!(adv_u,adv_v,q_u,q_v,V_u,U_v)

    #= Smagorinsky-like biharmonic diffusion
    Lu1 + Lu2 = dx[ νSmag dx(L(u))] + dy[ νSmag dy(L(u))]
    Lv1 + Lv2 = dx[ νSmag dx(L(v))] + dy[ νSmag dy(L(v))]
    =#
    @views νSmag_q .= (dudy .+ dvdx).^2     # reuse variable νSmag_q
    IqT!(shear,νSmag_q)
    @views νSmag .= cSmag*sqrt.((dudx .- dvdy).^2 .+ shear)
    ITq!(νSmag_q,νSmag)
    ∇²u!(Lu1,u)
    ∇²v!(Lv1,v)
    Gux!(dLudx,Lu1)
    Guy!(dLudy,Lu1)
    Gvx!(dLvdx,Lv1)
    Gvy!(dLvdy,Lv1)
    @views dLudx .= νSmag.*dLudx
    @views dLudy .= νSmag_q.*dLudy
    @views dLvdy .= νSmag.*dLvdy
    @views dLvdx .= νSmag_q.*dLvdx
    GTx!(Lu1,dLudx)        # reuse intermediate variable Lu1
    Gqy!(Lu2,dLudy)
    GTy!(Lv1,dLvdy)
    Gqx!(Lv2,dLvdx)

    # adding the terms
    momentum_u!(du,adv_u,dpdx,Lu1,Lu2,Fx)
    momentum_v!(dv.adv_v,dpdy,Lv1,Lv2)
    continuity!(dη,dUdx,dVdy)
end

function Uflux_nonperiodic!(U,u,h_u)
    @views U .= u[2:end-1,2:end-1].*h_u
end

function Uflux_periodic!(U,u,h_u)
    @views U .= u[3:end-1,2:end-1].*h_u
end

function Vflux!(V,v,h_v)
    @views V .= v[2:end-1,2:end-1].*h_v
end

function Bernoulli_nonperiodic!(p,KEu,KEv,η)
     @views p .= one_half*(KEu[:,2:end-1] .+ KEv[2:end-1,:]) .+ g*η
end

function Bernoulli_periodic!(p,KEu,KEv,η)
     @views p .= one_half*(KEu[2:end,2:end-1] .+ KEv[2:end-1,:]) .+ g*η
end


function PV_nonperiodic!(q,f_q,dvdx,dudy,h_q)
    #TODO this depends on the boundary conditions
    @views q .= (f_q .+ dvdx .- dudy) ./ h_q
end

function PV_periodic!(q,f_q,dvdx,dudy,h_q)
    #TODO this depends on the boundary conditions
    @views q .= (f_q .+ dvdx .- dudy) ./ h_q
end


function PV_adv!(adv_u,adv_v,q_u,q_v,V_u,U_v)
    @views adv_u .= q_u.*V_u
    @views adv_v .= -q_v.*U_v
end

function momentum_u!(du,adv_u,dpdx,Lu1,Lu2,Fx)
    @views du[3:end-2,3:end-2] .= adv_u .- dpdx[:,2:end-1] .+ Lu1 .+ Lu2 .+ Fx
end

function momentum_v!(dv,adv_v,dpdy,Lv1,Lv2)
    @views dv[3:end-2,3:end-2] .= adv_v .- dpdy[2:end-1,:] .+ Lv1.+Lv2
end

function continuity!(dη,dUdx,dVdy)
    # cut off the redundant halo and only copy the non-halo points into dη
    @inbounds for i ∈ 2:ny+1
        for j ∈ 2:nx+1
            dη[j,i] = -(dUdx[j-1,i] + dVdy[j,i-1])
        end
    end
end

if bc_x == "periodic"
    Uflux! = Uflux_periodic!
    Bernoulli! = Bernoulli_periodic!
else
    Uflux! = Uflux_nonperiodic!
    Bernoulli! = Bernoulli_nonperiodic!
end



# function continuity_mat!(dη,dUdx,dVdy)
#     # cut off the redundant halo and only copy the non-halo points into dη
#     @views @inbounds dη[2:end-1,2:end-1] .= -(dUdx[:,2:end-1] .+ dVdy[2:end-1,:])
# end

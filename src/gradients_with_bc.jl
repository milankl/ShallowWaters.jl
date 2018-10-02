# GRADIENT OPERATOR WITH DIMENSION 1/Δ

function GTx_nonperiodic!(dTdx::Matrix{Numtype},T::Matrix{Numtype})
    #= Calculates the gradient in x-direction on the T-grid.
    The result dTdx sits on the u-grid. =#

    dTdx[:,:] = one_over_Δ*(T[2:end,:] - T[1:end-1,:])
end

function GTx_periodic!(dTdx::Matrix{Numtype},T::Matrix{Numtype})
    #= Calculates the gradient in x-direction on the T-grid.
    The result dTdx sits on the u-grid. =#

    dTdx[1,:] = one_over_Δ*(T[1,:] - T[end,:])
    dTdx[2:end,:] = one_over_Δ*(T[2:end,:] - T[1:end-1,:])
end

function GTy!(dTdy::Matrix{Numtype},T::Matrix{Numtype})
    #= Calculates the gradient in y-direction on the T-grid.
    The result dTdy sits on the v-grid. =#

    dTdy[:,:] = one_over_Δ*(T[:,2:end] - T[:,1:end-1])
end

function Gux_nonperiodic!(dudx::Matrix{Numtype},u::Matrix{Numtype})
    #= Calculates the gradient in x-direction on the u-grid.
    The result dudx sits on the T-grid. =#

    dudx[1,:] = one_over_Δ*u[1,:]              # -0, kinematic bc
    dudx[2:end-1,:] = one_over_Δ*(u[2:end,:] - u[1:end-1,:])
    dudx[end,:] = (-one_over_Δ)*u[end,:]       # +0, kinematic bc
end

function Gux_periodic!(dudx::Matrix{Numtype},u::Matrix{Numtype})
    #= Calculates the gradient in x-direction on the u-grid.
    The result dudx sits on the T-grid. =#

    dudx[1:end-1,:] = one_over_Δ*(u[2:end,:]-u[1:end-1,:])
    dudx[end,:] = one_over_Δ*(u[1,:]-u[end,:])
end

function Gvy!(dvdy::Matrix{Numtype},v::Matrix{Numtype})
    #= Calculates the gradient in y-direction on the v-grid.
    The result dvdy sits on the T-grid. =#

    dvdy[:,2:end-1] = one_over_Δ*(v[:,2:end] - v[:,1:end-1])
    dvdy[:,1] = one_over_Δ*v[:,1]              # -0, kinematic bc
    dvdy[:,end] = (-one_over_Δ)*v[:,end]       # +0, kinematic bc
end

function Guy_nonperiodic!(dudy::Matrix{Numtype},u::Matrix{Numtype})
    #= Calculates the gradient in y-direction on the u-grid.
    The result dudy sits on the q-grid. =#

    dudy[2:end-1,2:end-1] = one_over_Δ*(u[:,2:end] - u[:,1:end-1])
    dudy[2:end-1,1] = α_over_Δ*u[:,1] #  α is the lateral boundary condition parameter
    dudy[2:end-1,end] = -α_over_Δ*u[:,end]
    dudy[1,:] .= zeero       #TODO probably redundant if initialised with zeros
    dudy[end,:] .= zeero
end

function Guy_periodic!(dudy::Matrix{Numtype},u::Matrix{Numtype})
    #= Calculates the gradient in y-direction on the u-grid.
    The result dudy sits on the q-grid. =#

    dudy[:,2:end-1] = one_over_Δ*(u[:,2:end] - u[:,1:end-1])
    dudy[:,1] = α_over_Δ*u[:,1] #  α is the lateral boundary condition parameter
    dudy[:,end] = -α_over_Δ*u[:,end]
end

function Gvx_nonperiodic!(dvdx::Matrix{Numtype},v::Matrix{Numtype})
    #= Calculates the gradient in x-direction on the v-grid.
    The result dvdx sits on the q-grid. =#

    dvdx[2:end-1,2:end-1] = one_over_Δ*(v[2:end,:] - v[1:end-1,:])
    dvdx[1,2:end-1] = α_over_Δ*v[1,:] #  α is the lateral boundary condition parameter
    dvdx[end,2:end-1] = -α_over_Δ*v[end,:]
    dvdx[:,1] .= zeero      #TODO redundant if initialised with zeros
    dvdx[:,end] .= zeero
end

function Gvx_periodic!(dvdx::Matrix{Numtype},v::Matrix{Numtype})
    #= Calculates the gradient in x-direction on the v-grid.
    The result dvdx sits on the q-grid. =#

    dvdx[2:end,2:end-1] = one_over_Δ*(v[2:end,:] - v[1:end-1,:])
    dvdx[1,2:end-1] = one_over_Δ*(v[1,:]-v[end,:])
    dvdx[:,1] .= zeero      #TODO redundant if initialised with zeros
    dvdx[:,end] .= zeero
end

function Gqx_nonperiodic!(dqdx::Matrix{Numtype},q::Matrix{Numtype})
    #= Calculates the gradient in x-direction on the q-grid.
    The result dqdx sits on the v-grid. =#

    dqdx[:,:] = one_over_Δ*(q[2:end,2:end-1] - q[1:end-1,2:end-1])
end

function Gqx_periodic!(dqdx::Matrix{Numtype},q::Matrix{Numtype})
    #= Calculates the gradient in x-direction on the q-grid.
    The result dqdx sits on the v-grid. =#

    dqdx[1:end-1,:] = one_over_Δ*(q[2:end,2:end-1] - q[1:end-1,2:end-1])
    dqdx[end,:] = one_over_Δ*(q[1,2:end-1] - q[end,2:end-1])
end

function Gqy_nonperiodic!(dqdy::Matrix{Numtype},q::Matrix{Numtype})
    #= Calculates the gradient in y-direction on the q-grid.
    The result dqdy sits on the u-grid. =#

    dqdy[:,:] = one_over_Δ*(q[2:end-1,2:end] - q[2:end-1,1:end-1])
end

function Gqy_periodic!(dqdy::Matrix{Numtype},q::Matrix{Numtype})
    #= Calculates the gradient in y-direction on the q-grid.
    The result dqdy sits on the u-grid. =#

    dqdy[:,:] = one_over_Δ*(q[:,2:end] - q[:,1:end-1])
end

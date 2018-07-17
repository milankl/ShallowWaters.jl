function ∇²u_nonperiodic!(du::AbstractMatrix,u::AbstractMatrix)
    #= ∇²u is the Laplace-operator d/dx^2 + d/dy^2 applied to a variable
    on the u-grid. The result du sits again on the u-grid.

    The 1/dx^2 factor is omitted and moved into the viscosity coefficient. =#

    # five-point stencil (1,1,-4,1,1) for the interior
    # central point + RIGHT + UP + LEFT + DOWN
    @views du[2:end-1,2:end-1] .= minus_4*u[2:end-1,2:end-1] .+ u[3:end,2:end-1] .+ u[2:end-1,3:end] .+ u[1:end-2,2:end-1] .+ u[2:end-1,1:end-2]

    # left edge, replace LEFT with zeros (kinematic BC)
    # central point .+ RIGHT + UP + DOWN
    @views du[1,2:end-1] .= minus_4*u[1,2:end-1] .+ u[2,2:end-1] .+ u[1,3:end] .+ u[1,1:end-2]

    # right edge, replace RIGHT with zeros (kinematic BC)
    # central point + LEFT + UP + DOWN
    @views du[end,2:end-1] .= minus_4*u[end,2:end-1] .+ u[end-1,2:end-1] .+ u[end,3:end] .+ u[end,1:end-2]

    # top edge (use boundary condition parameter α)
    # central point + LEFT + RIGHT + DOWN
    @views du[2:end-1,end] .= minus_3_minus_α*u[2:end-1,end] .+ u[1:end-2,end] .+ u[3:end,end] .+ u[2:end-1,end-1]

    # bottom edge (use boundary condition parameter α)
    # central point + DOWN + LEFT + RIGHT + UP
    @views du[2:end-1,1] .= minus_3_minus_α*u[2:end-1,1] .+ u[1:end-2,1] .+ u[3:end,1] .+ u[2:end-1,2]

    # bottom left (use boundary condition parameter α)
    @views du[1,1] .= minus_3_minus_α*u[1,1] .+ u[2,1] .+ u[1,2]

    # bottom right
    @views du[end,1] .= minus_3_minus_α*u[end,1] .+ u[end,2] .+ u[end-1,1]

    # top left
    @views du[1,end] .= minus_3_minus_α*u[1,end] .+ u[2,end] .+ u[1,end-1]

    # top right
    @views du[end,end] .= minus_3_minus_α*u[end,end] .+ u[end-1,end] .+ u[end,end-1]
end

function ∇²u_periodic!(du::AbstractMatrix,u::AbstractMatrix)
    #= ∇²u is the Laplace-operator d/dx^2 + d/dy^2 applied to a variable
    on the u-grid. The result du sits again on the u-grid.

    The 1/dx^2 factor is omitted and moved into the viscosity coefficient.=#

    # five-point stencil (1,1,-4,1,1) for the interior
    # central point + RIGHT + UP + LEFT + DOWN
    @views du[2:end-1,2:end-1] .= minus_4*u[2:end-1,2:end-1] .+ u[3:end,2:end-1] .+ u[2:end-1,3:end] .+ u[1:end-2,2:end-1] .+ u[2:end-1,1:end-2]

    # left edge - periodic BC
    # central point + RIGHT + UP + LEFT + DOWN
    @views du[1,2:end-1] .= minus_4*u[1,2:end-1] .+ u[2,2:end-1] .+ u[1,3:end] .+ u[end,2:end-1] .+ u[1,1:end-2]

    # right edge - periodic BC
    # central point + LEFT + UP + RIGHT + DOWN
    @views du[end,2:end-1] .= minus_4*u[end,2:end-1] .+ u[end-1,2:end-1] .+ u[end,3:end] .+ u[1,2:end-1] .+ u[end,1:end-2]

    # top edge (use boundary condition parameter α)
     # central point + TOP + LEFT + RIGHT + DOWN
    @views du[2:end-1,end] .= minus_3_minus_α*u[2:end-1,end] .+ u[1:end-2,end] .+ u[3:end,end] .+ u[2:end-1,end-1]

    # bottom edge (use boundary condition parameter α)
    # central point + DOWN + LEFT + RIGHT + UP
    @views du[2:end-1,1] .= minus_3_minus_α*u[2:end-1,1] .+ u[1:end-2,1] .+ u[3:end,1] .+ u[2:end-1,2]

    # bottom left
    @views du[1,1] .= minus_3_minus_α*u[1,1] .+ u[end,1] .+ u[2,1] .+ u[1,2]

    # # bottom right
    @views du[end,1] .= minus_3_minus_α*u[end,1] .+ u[1,1] .+ u[end,2] .+ u[end-1,1]

    # top left
    @views du[1,end] .= minus_3_minus_α*u[1,end] .+ u[end,end] .+ u[2,end] .+ u[1,end-1]

    # top right
    @views du[end,end] .= minus_3_minus_α*u[end,end] .+ u[1,end] .+ u[end-1,end] .+ u[end,end-1]
end

function ∇²v_nonperiodic!(dv::AbstractMatrix,v::AbstractMatrix)
    #= ∇²v is the Laplace-operator d/dx^2 + d/dy^2 applied to a variable
    on the v-grid. The result dv sits again on the v-grid.

    The 1/dx^2 factor is omitted and moved into the viscosity coefficient.=#

    # five-point stencil (1,1,-4,1,1) for the interior
    # central point + RIGHT + UP + LEFT + DOWN
    @views dv[2:end-1,2:end-1] .= minus_4*v[2:end-1,2:end-1] .+ v[3:end,2:end-1] .+ v[2:end-1,3:end] .+ v[1:end-2,2:end-1] .+ v[2:end-1,1:end-2]

    # left edge, use boundary condition parameter α
    # central point + RIGHT + UP + DOWN
    @views dv[1,2:end-1] .= minus_3_minus_α*v[1,2:end-1] .+ v[2,2:end-1] .+ v[1,3:end] .+ v[1,1:end-2]

    # right edge, use boundary condition parameter α
    # central point + LEFT + UP  + DOWN
    @views dv[end,2:end-1] .= minus_3_minus_α*v[end,2:end-1] .+ v[end-1,2:end-1] .+ v[end,3:end] .+ v[end,1:end-2]

    # top edge (use kinematic BC)
    # central point + TOP  # LEFT + RIGHT # DOWN
    @views dv[2:end-1,end] .= minus_4*v[2:end-1,end] .+ v[1:end-2,end] .+ v[3:end,end] .+ v[2:end-1,end-1]

    # bottom edge (use kinematic BC)
    # central point + DOWN  # LEFT + RIGHT # UP
    @views dv[2:end-1,1] .= minus_4*v[2:end-1,1] .+ v[1:end-2,1] .+ v[3:end,1] .+ v[2:end-1,2]

    # bottom left (use boundary condition parameter α)
    @views dv[1,1] .= minus_3_minus_α*v[1,1] .+ v[2,1]  .+ v[1,2]

    # bottom right
    @views dv[end,1] .= minus_3_minus_α*v[end,1] .+ v[end,2] .+ v[end-1,1]

    # top left
    @views dv[1,end] .= minus_3_minus_α*v[1,end] .+ v[2,end] .+ v[1,end-1]

    # top right
    @views dv[end,end] .= minus_3_minus_α*v[end,end] .+ v[end-1,end] .+ v[end,end-1]
end

function ∇²v_periodic!(dv::AbstractMatrix,v::AbstractMatrix)
    #= ∇²v is the Laplace-operator d/dx^2 + d/dy^2 applied to a variable
    on the v-grid. The result dv sits again on the v-grid.

    The 1/dx^2 factor is omitted and moved into the viscosity coefficient.=#

    # five-point stencil (1,1,-4,1,1) for the interior
    # central point # RIGHT + UP # LEFT + DOWN
    @views dv[2:end-1,2:end-1] .= minus_4*v[2:end-1,2:end-1] .+ v[3:end,2:end-1] .+ v[2:end-1,3:end] .+ v[1:end-2,2:end-1] .+ v[2:end-1,1:end-2]

    # left edge - periodic BC
    # central point # RIGHT + UP  # LEFT + DOWN
    @views dv[1,2:end-1] .= minus_4*v[1,2:end-1] .+ v[2,2:end-1] .+ v[1,3:end] .+ v[end,2:end-1] .+ v[1,1:end-2]

    # right edge - periodic BC
    # central point + LEFT + UP + RIGHT + DOWN
    @views dv[end,2:end-1] .= minus_4*v[end,2:end-1] .+ v[end-1,2:end-1] .+ v[end,3:end] .+  v[1,2:end-1] .+ v[end,1:end-2]

    # top edge (use kinematic BC)
    # central point + TOP  # LEFT + RIGHT # DOWN
    @views dv[2:end-1,end] .= minus_4*v[2:end-1,end] .+ v[1:end-2,end] .+ v[3:end,end] .+ v[2:end-1,end-1]

    # bottom edge (use kinematic BC)
    # central point + DOWN + LEFT + RIGHT + UP
    @views dv[2:end-1,1] .= minus_4*v[2:end-1,1] .+ v[1:end-2,1] .+ v[3:end,1] .+ v[2:end-1,2]

    # bottom left
    @views dv[1,1] .= minus_3_minus_α*v[1,1] .+ v[end,1] .+ v[2,1]  .+ v[1,2]

    # bottom right
    @views dv[end,1] .= minus_3_minus_α*v[end,1] .+ v[1,1] .+ v[end,2] .+ v[end-1,1]

    # top left
    @views dv[1,end] .= minus_3_minus_α*v[1,end] .+ v[end,end] .+ v[2,end] .+ v[1,end-1]

    # top right
    @views dv[end,end] .= minus_3_minus_α*v[end,end] .+ v[1,end] .+ v[end-1,end] .+ v[end,end-1]
end

function ∇⁴u_periodic!(du::AbstractMatrix,du2::AbstractMatrix,u::AbstractMatrix)
    #= ∇⁴u is the biharmonic operator d/dx^4 + d/dy^4 + 2*(d/dx^2)*(d/dy^2)
    applied to a variable on the u-grid. The result du sits again on the u-grid. =#

    # apply the harmonic Laplacian twice
    ∇²u_periodic!(du,u)
    ∇²u_periodic!(du2,du)
end

function ∇⁴u_nonperiodic!(du::AbstractMatrix,du2::AbstractMatrix,u::AbstractMatrix)
    #= ∇⁴u is the biharmonic operator d/dx^4 + d/dy^4 + 2*(d/dx^2)*(d/dy^2)
    applied to a variable on the u-grid. The result du sits again on the u-grid. =#

    # apply the harmonic Laplacian twice
    ∇²u_nonperiodic!(du,u)
    ∇²u_nonperiodic!(du2,du)
end

function ∇⁴v_periodic!(dv::AbstractMatrix,dv2::AbstractMatrix,v::AbstractMatrix)
    #= ∇⁴v is the biharmonic operator d/dx^4 + d/dy^4 + 2*(d/dx^2)*(d/dy^2)
    applied to a variable on the v-grid. The result dv sits again on the v-grid. =#

    # apply the harmonic Laplacian twice
    ∇²v_periodic!(dv,v)
    ∇²v_periodic!(dv2,v)
end

function ∇⁴v_nonperiodic!(dv::AbstractMatrix,dv2::AbstractMatrix,v::AbstractMatrix)
    #= ∇⁴v is the biharmonic operator d/dx^4 + d/dy^4 + 2*(d/dx^2)*(d/dy^2)
    applied to a variable on the v-grid. The result dv sits again on the v-grid. =#

    # apply the harmonic Laplacian twice
    ∇²v_nonperiodic!(dv,v)
    ∇²v_nonperiodic!(dv2,v)
end

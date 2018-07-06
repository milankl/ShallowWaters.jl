function ∇²u_nonperiodic!(du::Matrix{Numtype},u::Matrix{Numtype})
    #= ∇²u is the Laplace-operator d/dx^2 + d/dy^2 applied to a variable
    on the u-grid. The result du sits again on the u-grid.

    The 1/dx^2 factor is omitted and moved into the viscosity coefficient. =#

    # five-point stencil (1,1,-4,1,1) for the interior
    du[2:end-1,2:end-1] = minus_4*u[2:end-1,2:end-1] +             # central point
                          u[3:end,2:end-1] + u[2:end-1,3:end] +       # RIGHT + UP
                          u[1:end-2,2:end-1] + u[2:end-1,1:end-2]     # LEFT + DOWN

    # left edge, replace LEFT with zeros (kinematic BC)
    du[1,2:end-1] = minus_4*u[1,2:end-1] +                         # central point
                    u[2,2:end-1] + u[1,3:end] +                       # RIGHT + UP
                    u[1,1:end-2]                                      # DOWN

    # right edge, replace RIGHT with zeros (kinematic BC)
    du[end,2:end-1] = minus_4*u[end,2:end-1] +                     # central point
                      u[end-1,2:end-1] + u[end,3:end] +               # LEFT + UP
                      u[end,1:end-2]                                  # DOWN

    # top edge (use boundary condition parameter α)
    du[2:end-1,end] = minus_3_minus_α*u[2:end-1,end] +            # central point + TOP
                      u[1:end-2,end] + u[3:end,end] +                 # LEFT + RIGHT
                      u[2:end-1,end-1]                                # DOWN

    # bottom edge (use boundary condition parameter α)
    du[2:end-1,1] = minus_3_minus_α*u[2:end-1,1] +                # central point + DOWN
                    u[1:end-2,1] + u[3:end,1] +                       # LEFT + RIGHT
                    u[2:end-1,2]                                      # UP

    # bottom left (use boundary condition parameter α)
    du[1,1] = minus_3_minus_α*u[1,1] + u[2,1] + u[1,2]

    # bottom right
    du[end,1] = minus_3_minus_α*u[end,1] + u[end,2] + u[end-1,1]

    # top left
    du[1,end] = minus_3_minus_α*u[1,end] + u[2,end] + u[1,end-1]

    # top right
    du[end,end] = minus_3_minus_α*u[end,end]+u[end-1,end] + u[end,end-1]
end

function ∇²u_periodic!(du::Matrix{Numtype},u::Matrix{Numtype})
    #= ∇²u is the Laplace-operator d/dx^2 + d/dy^2 applied to a variable
    on the u-grid. The result du sits again on the u-grid.

    The 1/dx^2 factor is omitted and moved into the viscosity coefficient.=#

    # five-point stencil (1,1,-4,1,1) for the interior
    du[2:end-1,2:end-1] = minus_4*u[2:end-1,2:end-1] +               # central point
                            u[3:end,2:end-1] + u[2:end-1,3:end] +       # RIGHT + UP
                            u[1:end-2,2:end-1] + u[2:end-1,1:end-2]     # LEFT + DOWN

    # left edge - periodic BC
    du[1,2:end-1] = minus_4*u[1,2:end-1] +                           # central point
                            u[2,2:end-1] + u[1,3:end] +                 # RIGHT + UP
                            u[end,2:end-1] + u[1,1:end-2]               # LEFT + DOWN

    # right edge - periodic BC
    du[end,2:end-1] = minus_4*u[end,2:end-1] +                       # central point
                            u[end-1,2:end-1] + u[end,3:end] +           # LEFT + UP
                            u[1,2:end-1] + u[end,1:end-2]               # RIGHT + DOWN

    # top edge (use boundary condition parameter α)
    du[2:end-1,end] = minus_3_minus_α*u[2:end-1,end] +              # central point + TOP
                            u[1:end-2,end] + u[3:end,end] +             # LEFT + RIGHT
                            u[2:end-1,end-1]                            # DOWN

    # bottom edge (use boundary condition parameter α)
    du[2:end-1,1] = minus_3_minus_α*u[2:end-1,1] +                  # central point + DOWN
                            u[1:end-2,1] + u[3:end,1] +                 # LEFT + RIGHT
                            u[2:end-1,2]                                # UP

    # bottom left
    du[1,1] = minus_3_minus_α*u[1,1] + u[end,1] + u[2,1] + u[1,2]

    # # bottom right
    du[end,1] = minus_3_minus_α*u[end,1] + u[1,1] + u[end,2] + u[end-1,1]

    # top left
    du[1,end] = minus_3_minus_α*u[1,end] + u[end,end] + u[2,end] + u[1,end-1]

    # top right
    du[end,end] = minus_3_minus_α*u[end,end] + u[1,end] + u[end-1,end] + u[end,end-1]
end

function ∇²v_nonperiodic!(dv::Matrix{Numtype},v::Matrix{Numtype})
    #= ∇²v is the Laplace-operator d/dx^2 + d/dy^2 applied to a variable
    on the v-grid. The result dv sits again on the v-grid.

    The 1/dx^2 factor is omitted and moved into the viscosity coefficient.=#

    # five-point stencil (1,1,-4,1,1) for the interior
    dv[2:end-1,2:end-1] = minus_4*v[2:end-1,2:end-1] +               # central point
                            v[3:end,2:end-1] + v[2:end-1,3:end] +       # RIGHT + UP
                            v[1:end-2,2:end-1] + v[2:end-1,1:end-2]     # LEFT + DOWN

    # left edge, use boundary condition parameter α
    dv[1,2:end-1] = minus_3_minus_α*v[1,2:end-1] +                  # central point
                            v[2,2:end-1] + v[1,3:end] +                 # RIGHT + UP
                            v[1,1:end-2]                                # DOWN

    # right edge, use boundary condition parameter α
    dv[end,2:end-1] = minus_3_minus_α*v[end,2:end-1] +              # central point
                            v[end-1,2:end-1] + v[end,3:end] +           # LEFT + UP
                            v[end,1:end-2]                              # DOWN

    # top edge (use kinematic BC)
    dv[2:end-1,end] = minus_4*v[2:end-1,end] +                       # central point + TOP
                            v[1:end-2,end] + v[3:end,end] +             # LEFT + RIGHT
                            v[2:end-1,end-1]                            # DOWN

    # bottom edge (use kinematic BC)
    dv[2:end-1,1] = minus_4*v[2:end-1,1] +                           # central point + DOWN
                            v[1:end-2,1] + v[3:end,1] +                 # LEFT + RIGHT
                            v[2:end-1,2]                                # UP

    # bottom left (use boundary condition parameter α)
    dv[1,1] = minus_3_minus_α*v[1,1] + v[2,1]  + v[1,2]

    # bottom right
    dv[end,1] = minus_3_minus_α*v[end,1] + v[end,2] + v[end-1,1]

    # top left
    dv[1,end] = minus_3_minus_α*v[1,end] + v[2,end] + v[1,end-1]

    # top right
    dv[end,end] = minus_3_minus_α*v[end,end] + v[end-1,end] + v[end,end-1]
end

function ∇²v_periodic!(dv::Matrix{Numtype},v::Matrix{Numtype})
    #= ∇²v is the Laplace-operator d/dx^2 + d/dy^2 applied to a variable
    on the v-grid. The result dv sits again on the v-grid.

    The 1/dx^2 factor is omitted and moved into the viscosity coefficient.=#

    # five-point stencil (1,1,-4,1,1) for the interior
    dv[2:end-1,2:end-1] = minus_4*v[2:end-1,2:end-1] +               # central point
                            v[3:end,2:end-1] + v[2:end-1,3:end] +       # RIGHT + UP
                            v[1:end-2,2:end-1] + v[2:end-1,1:end-2]     # LEFT + DOWN

    # left edge - periodic BC
    dv[1,2:end-1] = minus_4*v[1,2:end-1] +                           # central point
                            v[2,2:end-1] + v[1,3:end] +                 # RIGHT + UP
                            v[end,2:end-1] + v[1,1:end-2]               # LEFT + DOWN

    # right edge - periodic BC
    dv[end,2:end-1] = minus_4*v[end,2:end-1] +                       # central point
                            v[end-1,2:end-1] + v[end,3:end] +           # LEFT + UP
                            v[1,2:end-1] + v[end,1:end-2]               # RIGHT + DOWN

    # top edge (use kinematic BC)
    dv[2:end-1,end] = minus_4*v[2:end-1,end] +                       # central point + TOP
                            v[1:end-2,end] + v[3:end,end] +             # LEFT + RIGHT
                            v[2:end-1,end-1]                            # DOWN

    # bottom edge (use kinematic BC)
    dv[2:end-1,1] = minus_4*v[2:end-1,1] +                           # central point + DOWN
                            v[1:end-2,1] + v[3:end,1] +                 # LEFT + RIGHT
                            v[2:end-1,2]                                # UP

    # bottom left
    dv[1,1] = minus_3_minus_α*v[1,1] + v[end,1] + v[2,1]  + v[1,2]

    # bottom right
    dv[end,1] = minus_3_minus_α*v[end,1] + v[1,1] + v[end,2] + v[end-1,1]

    # top left
    dv[1,end] = minus_3_minus_α*v[1,end] + v[end,end] + v[2,end] + v[1,end-1]

    # top right
    dv[end,end] = minus_3_minus_α*v[end,end] + v[1,end] + v[end-1,end] + v[end,end-1]
end

function ∇⁴u_periodic!(du::Matrix{Numtype},du2::Matrix{Numtype},u::Matrix{Numtype})
    #= ∇⁴u is the biharmonic operator d/dx^4 + d/dy^4 + 2*(d/dx^2)*(d/dy^2)
    applied to a variable on the u-grid. The result du sits again on the u-grid. =#

    # apply the harmonic Laplacian twice
    ∇²u_periodic!(du,u)
    ∇²u_periodic!(du2,du)
end

function ∇⁴u_nonperiodic!(du::Matrix{Numtype},du2::Matrix{Numtype},u::Matrix{Numtype})
    #= ∇⁴u is the biharmonic operator d/dx^4 + d/dy^4 + 2*(d/dx^2)*(d/dy^2)
    applied to a variable on the u-grid. The result du sits again on the u-grid. =#

    # apply the harmonic Laplacian twice
    ∇²u_nonperiodic!(du,u)
    ∇²u_nonperiodic!(du2,du)
end

function ∇⁴v_periodic!(dv::Matrix{Numtype},dv2::Matrix{Numtype},v::Matrix{Numtype})
    #= ∇⁴v is the biharmonic operator d/dx^4 + d/dy^4 + 2*(d/dx^2)*(d/dy^2)
    applied to a variable on the v-grid. The result dv sits again on the v-grid. =#

    # apply the harmonic Laplacian twice
    ∇²v_periodic!(dv,v)
    ∇²v_periodic!(dv2,v)
end

function ∇⁴v_nonperiodic!(dv::Matrix{Numtype},dv2::Matrix{Numtype},v::Matrix{Numtype})
    #= ∇⁴v is the biharmonic operator d/dx^4 + d/dy^4 + 2*(d/dx^2)*(d/dy^2)
    applied to a variable on the v-grid. The result dv sits again on the v-grid. =#

    # apply the harmonic Laplacian twice
    ∇²v_nonperiodic!(dv,v)
    ∇²v_nonperiodic!(dv2,v)
end

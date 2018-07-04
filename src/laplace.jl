function Lu_nonperiodic(du,u)
    #= Lu is the Laplace-operator d/dx^2 + d/dy^2 applied to a variable
    on the u-grid. The result du sits again on the u-grid.

    The 1/dx^2 factor is omitted and moved into the viscosity coefficient. =#

    # five-point stencil (1,1,-4,1,1) for the interior
    du[2:end-1,2:end-1] = -four*u[2:end-1,2:end-1] +                  # central point
                          u[3:end,2:end-1] + u[2:end-1,3:end] +       # RIGHT + UP
                          u[1:end-2,2:end-1] + u[2:end-1,1:end-2]     # LEFT + DOWN

    # left edge, replace LEFT with zeros (kinematic BC)
    du[1,2:end-1] = -four*u[1,2:end-1] +                              # central point
                    u[2,2:end-1] + u[1,3:end] +                       # RIGHT + UP
                    u[1,1:end-2]                                      # DOWN

    # right edge, replace RIGHT with zeros (kinematic BC)
    du[end,2:end-1] = -four*u[end,2:end-1] +                          # central point
                      u[end-1,2:end-1] + u[end,3:end] +               # LEFT + UP
                      u[end,1:end-2]                                  # DOWN

    # top edge (use boundary condition parameter α)
    du[2:end-1,end] = -four_plus_α*u[2:end-1,end] +                   # central point + TOP
                      u[1:end-2,end] + u[3:end,end] +                 # LEFT + RIGHT
                      u[2:end-1,end-1]                                # DOWN

    # bottom edge (use boundary condition parameter α)
    du[2:end-1,1] = -four_minus_α*u[2:end-1,1] +                      # central point + DOWN
                    u[1:end-2,1] + u[3:end,1] +                       # LEFT + RIGHT
                    u[2:end-1,2]                                      # UP

    # bottom left (use boundary condition parameter α)
    du[1,1] = -four_minus_α*u[1,1] + u[2,1] + u[1,2]                  # central + RIGHT + UP

    # bottom right
    du[end,1] = -four_minus_α*u[end,1] + u[end,2] + u[end-1,1]        # central + UP + LEFT

    # top left
    du[1,end] = -four_plus_α*u[1,end] + u[2,end] + u[1,end-1]         # central + RIGHT + DOWN

    # top right
    du[end,end] = -four_plus_α*u[end,end]+u[end-1,end] + u[end,end-1] # central + LEFT + DOWN

    return du
end

function Lu_periodic(du,u)
    #= Lu is the Laplace-operator d/dx^2 + d/dy^2 applied to a variable
    on the u-grid. The result du sits again on the u-grid.

    The 1/dx^2 factor is omitted and moved into the viscosity coefficient.=#

    # five-point stencil (1,1,-4,1,1) for the interior
    du[2:end-1,2:end-1] = -four*u[2:end-1,2:end-1] +      # central point
                            u[3:end,2:end-1] + u[2:end-1,3:end] +       # RIGHT + UP
                            u[1:end-2,2:end-1] + u[2:end-1,1:end-2]    # LEFT + DOWN

    # left edge - periodic BC
    du[1,2:end-1] = -four*u[1,2:end-1] +                  # central point
                            u[2,2:end-1] + u[1,3:end] +                 # RIGHT + UP
                            u[end,2:end-1] + u[1,1:end-2]              # LEFT + DOWN

    # right edge - periodic BC
    du[end,2:end-1] = -four*u[end,2:end-1] +              # central point
                            u[end-1,2:end-1] + u[end,3:end] +           # LEFT + UP
                            u[1,2:end-1] + u[end,1:end-2]              # RIGHT + DOWN

    # top edge (use boundary condition parameter α)
    du[2:end-1,end] = -four_plus_α*u[2:end-1,end] +   # central point + TOP
                            u[1:end-2,end] + u[3:end,end] +        # LEFT + RIGHT
                            u[2:end-1,end-1]                      # DOWN

    # bottom edge (use boundary condition parameter α)
    du[2:end-1,1] = -four_minus_α*u[2:end-1,1] +   # central point + DOWN
                            u[1:end-2,1] + u[3:end,1] +        # LEFT + RIGHT
                            u[2:end-1,2]                      # UP

    # bottom left
    du[1,1] = -four_minus_α*u[1,1] +
                            u[end,1] + u[2,1]  + u[1,2]       # LEFT + RIGHT + UP
    # bottom right
    du[end,1] = -four_minus_α*u[end,1] +
                            u[1,1] + u[end,2] + u[end-1,1]    # RIGHT + UP + LEFT

    # top left
    du[1,end] = -four_plus_α*u[1,end] +
                            u[end,end] + u[2,end] + u[1,end-1]  # LEFT + RIGHT + DOWN

    # top right
    du[end,end] = -four_plus_α*u[end,end] +
                            u[1,end] + u[end-1,end] + u[end,end-1] # RIGHT + LEFT + DOWN
    return du
end

function Lv_nonperiodic(dv,v)
    #= Lv is the Laplace-operator d/dx^2 + d/dy^2 applied to a variable
    on the v-grid. The result dv sits again on the v-grid.

    The 1/dx^2 factor is omitted and moved into the viscosity coefficient.=#

    # five-point stencil (1,1,-4,1,1) for the interior
    dv[2:end-1,2:end-1] = -four*v[2:end-1,2:end-1] +      # central point
                            v[3:end,2:end-1] + v[2:end-1,3:end] +       # RIGHT + UP
                            v[1:end-2,2:end-1] + v[2:end-1,1:end-2]    # LEFT + DOWN

    # left edge, use boundary condition parameter α
    dv[1,2:end-1] = -four_minus_α*v[1,2:end-1] +            # central point
                            v[2,2:end-1] + v[1,3:end] +                 # RIGHT + UP
                            v[1,1:end-2]                               # DOWN

    # right edge, use boundary condition parameter α
    dv[end,2:end-1] = -four_plus_α*v[end,2:end-1] +        # central point
                            v[end-1,2:end-1] + v[end,3:end] +           # LEFT + UP
                            v[end,1:end-2]                             # DOWN

    # top edge (use kinematic BC)
    dv[2:end-1,end] = -four*v[2:end-1,end] +         # central point + TOP
                            v[1:end-2,end] + v[3:end,end] +        # LEFT + RIGHT
                            v[2:end-1,end-1]                      # DOWN

    # bottom edge (use kinematic BC)
    dv[2:end-1,1] = -four*v[2:end-1,1] +         # central point + DOWN
                            v[1:end-2,1] + v[3:end,1] +        # LEFT + RIGHT
                            v[2:end-1,2]                      # UP

    # bottom left (use boundary condition parameter α)
    dv[1,1] = -four_minus_α*v[1,1] +
                            v[2,1]  + v[1,2]       # RIGHT + UP
    # bottom right
    dv[end,1] = -four_minus_α*v[end,1] +
                            v[end,2] + v[end-1,1]    # UP + LEFT

    # top left
    dv[1,end] = -four_plus_α*v[1,end] +
                            v[2,end] + v[1,end-1]  # RIGHT + DOWN

    # top right
    dv[end,end] = -four_plus_α*v[end,end] +
                            v[end-1,end] + v[end,end-1] # LEFT + DOWN
    return dv
end

function Lv_periodic(dv,v)
    #= Lv is the Laplace-operator d/dx^2 + d/dy^2 applied to a variable
    on the v-grid. The result dv sits again on the v-grid.

    The 1/dx^2 factor is omitted and moved into the viscosity coefficient.=#

    # five-point stencil (1,1,-4,1,1) for the interior
    dv[2:end-1,2:end-1] = -four*v[2:end-1,2:end-1] +      # central point
                            v[3:end,2:end-1] + v[2:end-1,3:end] +       # RIGHT + UP
                            v[1:end-2,2:end-1] + v[2:end-1,1:end-2]    # LEFT + DOWN

    # left edge - periodic BC
    dv[1,2:end-1] = -four*v[1,2:end-1] +                  # central point
                            v[2,2:end-1] + v[1,3:end] +                 # RIGHT + UP
                            v[end,2:end-1] + v[1,1:end-2]              # LEFT + DOWN

    # right edge - periodic BC
    dv[end,2:end-1] = -four*v[end,2:end-1] +              # central point
                            v[end-1,2:end-1] + v[end,3:end] +           # LEFT + UP
                            v[1,2:end-1] + v[end,1:end-2]              # RIGHT + DOWN

    # top edge (use kinematic BC)
    dv[2:end-1,end] = -four*v[2:end-1,end] +         # central point + TOP
                            v[1:end-2,end] + v[3:end,end] +        # LEFT + RIGHT
                            v[2:end-1,end-1]                      # DOWN

    # bottom edge (use kinematic BC)
    dv[2:end-1,1] = -four*v[2:end-1,1] +         # central point + DOWN
                            v[1:end-2,1] + v[3:end,1] +        # LEFT + RIGHT
                            v[2:end-1,2]                      # UP

    # bottom left
    dv[1,1] = -four_minus_α*v[1,1] +
                            v[end,1] + v[2,1]  + v[1,2]       # LEFT + RIGHT + UP
    # bottom right
    dv[end,1] = -four_minus_α*v[end,1] +
                            v[1,1] + v[end,2] + v[end-1,1]    # RIGHT + UP + LEFT

    # top left
    dv[1,end] = -four_plus_α*v[1,end] +
                            v[end,end] + v[2,end] + v[1,end-1]  # LEFT + RIGHT + DOWN

    # top right
    dv[end,end] = -four_plus_α*v[end,end] +
                            v[1,end] + v[end-1,end] + v[end,end-1] # RIGHT + LEFT + DOWN
    return dv
end

function LLu_periodic(du,du2,u)
    #= LLu is the biharmonic operator d/dx^4 + d/dy^4 + 2*(d/dx^2)*(d/dy^2)
    applied to a variable on the u-grid. The result du sits again on the u-grid. =#

    # apply the harmonic Laplacian twice
    du[:] = Lu_periodic(du,u)
    du2[:] = Lu_periodic(du2,du)

    return du2
end

function LLu_nonperiodic(du,du2,u)
    #= LLu is the biharmonic operator d/dx^4 + d/dy^4 + 2*(d/dx^2)*(d/dy^2)
    applied to a variable on the u-grid. The result du sits again on the u-grid. =#

    # apply the harmonic Laplacian twice
    du[:] = Lu_nonperiodic(du,u)
    du2[:] = Lu_nonperiodic(du2,du)

    return du2
end

function LLv_periodic(dv,dv2,v)
    #= LLv is the biharmonic operator d/dx^4 + d/dy^4 + 2*(d/dx^2)*(d/dy^2)
    applied to a variable on the v-grid. The result dv sits again on the v-grid. =#

    # apply the harmonic Laplacian twice
    dv[:] = Lv_periodic(dv,v)
    dv2[:] = Lv_periodic(dv2,v)

    return dv2
end

function LLv_nonperiodic(dv,dv2,v)
    #= LLv is the biharmonic operator d/dx^4 + d/dy^4 + 2*(d/dx^2)*(d/dy^2)
    applied to a variable on the v-grid. The result dv sits again on the v-grid. =#

    # apply the harmonic Laplacian twice
    dv[:] = Lv_nonperiodic(dv,v)
    dv2[:] = Lv_nonperiodic(dv2,v)

    return dv2
end

# renames the operator functions according to the boundary conditions set
# by bc_x

if bc_x == "periodic"    # i.e. periodic

    # GRADIENTS
    GTx! = GTx_periodic!
    # GTy same for periodic/nonperiodic
    Gux! = Gux_periodic!
    # Gvy same for periodic/nonperiodic
    Guy! = Guy_periodic!
    Gvx! = Gvx_periodic!
    Gqx! = Gqx_periodic!
    Gqy! = Gqy_periodic!

    # INTERPOLATIONS
    ITu! = ITu_periodic!
    # ITv same for periodic/nonperiodic
    IuT! = IuT_periodic!
    # IvT same for periodic/nonperiodic
    Iqu! = Iqu_periodic!
    Iqv! = Iqv_periodic!
    Ivq! = Ivq_periodic!
    Iuq! = Iuq_periodic!
    ITq! = ITq_periodic!
    IqT! = IqT_periodic!
    Iuv! = Iuv_periodic!
    Ivu! = Ivu_periodic!

    # Laplacians - always dimensionless
    ∇²u! = ∇²u_periodic!
    ∇²v! = ∇²v_periodic!
    ∇⁴u! = ∇⁴u_periodic!
    ∇⁴v! = ∇⁴v_periodic!

else            # i.e. non-periodic

    # GRADIENTS
    GTx! = GTx_nonperiodic!
    # GTy same for periodic/nonperiodic
    Gux! = Gux_nonperiodic!
    # Gvy same for periodic/nonperiodic
    Guy! = Guy_nonperiodic!
    Gvx! = Gvx_nonperiodic!
    Gqx! = Gqx_nonperiodic!
    Gqy! = Gqy_nonperiodic!

    # INTERPOLATIONS
    ITu! = ITu_nonperiodic!
    # ITv same for periodic/nonperiodic
    IuT! = IuT_nonperiodic!
    # IvT same for periodic/nonperiodic
    Iqu! = Iqu_nonperiodic!
    Iqv! = Iqv_nonperiodic!
    Ivq! = Ivq_nonperiodic!
    Iuq! = Iuq_nonperiodic!
    ITq! = ITq_nonperiodic!
    IqT! = IqT_nonperiodic!
    Iuv! = Iuv_nonperiodic!
    Ivu! = Ivu_nonperiodic!

    # Laplacians - always dimensionless
    ∇²u! = ∇²u_nonperiodic!
    ∇²v! = ∇²v_nonperiodic!
    ∇⁴u! = ∇⁴u_nonperiodic!
    ∇⁴v! = ∇⁴v_nonperiodic!
end

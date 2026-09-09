"""
This is an empty function to serve as a placeholder when Lux is not included and thus the extension is not included. 
"""

function CNN_momentum end

"""
The actual CNN_momentum function, commented out below, is in the extension ShallowWatersNNForcing.jl. To activate it and not
the empty function, one runs

    using ShallowWaters
    using Lux

which should setup the model to use the extension, and thus the actual NN forcing in the integration. From testing it's super 
important to Julia that you run 'using ShallowWaters' before running 'using Lux'. If they're run in the opposite order, then
the model structure CNNVars will remain empty, and the actual CNN_momentum function from the extension won't get run. However,
even with this just note that Julia's extensions are mildly finnicky (at least to me).
"""

# The function CNN_momentum takes 
#     u: the x-direction velocity
#     v: the y-direction velocity
#     S: the model defined in ShallowWaters
# The function returns the forcing tensor S, built from spatial derivatives of the 
# tensors
#     T_11,
#     T_12,
# and
#     T_22
# T_11, T_12, and T_22 are output from two separate CNNs
# The inputs to the CNNs will be the vorticity, shear, and stretch deformation fields,
# following the Zanna & Bolton parameterization contained in zanna_bolton_forcing.

# As of 09/08/26 the weights in this forcing function have been successfully trained in a flat bottom,
# barotropic gyre setup with the following parameters:
# ShallowWaters.Parameter(T=T;
#     output=false,
#     L_ratio=1,
#     g=9.81,
#     H=500,
#     wind_forcing_x="double_gyre",
#     Lx=3840e3,
#     seasonal_wind_x=false,
#     topography="flat",
#     bc="nonperiodic",
#     bottom_drag="quadratic",
#     tracer_advection=false,
#     tracer_relaxation=false,
#     nn_forcing_dissipation=true,
#     N=1,
#     α=2,
#     nx=128
# )
# This is not to say that the CNNs can't be trained in a different setup, but this is all that's been tested.
# Also, just like with the Zanna-Bolton parameterization, CNN_forcing_momentum has issues on a periodic domain.

# Warning about this function: the default weights that Lux chooses *will* cause the model to diverge, one needs to use
# tuned weights/tune the weights themselves to make the forcing work stably alongside model integration
# """
# function CNN_momentum(u, v, S)

#     Diag = S.Diag

#     # I think this needs to stay a Float32 because everything in Lux is Float32,
#     # if I remember correct there's an error that occurs without this line
#     T = Float32

#     @unpack nqx, nqy = Diag.CNNVars
#     @unpack dudx, dudy, dvdx, dvdy = Diag.CNNVars
#     @unpack γ₀, ζ, D, Dhat, Dhatq = Diag.CNNVars
#     @unpack model_Su, model_Sv = Diag.CNNVars
#     @unpack Su_layers, Sv_layers = Diag.CNNVars
#     @unpack Dhatq, ζT, DT, DhatT = Diag.CNNVars

#     @unpack T11, T22, T12 = Diag.CNNVars
#     @unpack dT11dx, dT12dy, dT12dx, dT22dy = Diag.CNNVars

#     @unpack res_Su, res_Sv = Diag.CNNVars
#     @unpack S_u, S_v = Diag.CNNVars

#     @unpack Δ, scale, f₀ = S.grid
#     @unpack halo, haloη, ep, nux, nuy, nvx, nvy = S.grid

#     mq,nq = size(ζ)
#     mTh,nTh = size(Dhat)
#     nx, ny = size(DT)

#     κ_BT = - γ₀ * Δ^2

#     ∂x!(dudx, u)
#     ∂y!(dudy, u)

#     ∂x!(dvdx, v)
#     ∂y!(dvdy, v)

#     # Relative vorticity and shear deformation, cell corners
#     @inbounds for j ∈ 1:nq
#         for k ∈ 1:mq
#             ζ[k,j] = dvdx[k+1,j+1] - dudy[k+1,j+1]
#             D[k,j] = dudy[k+1,j+1] + dvdx[k+1,j+1]
#         end
#     end

#     # Stretch deformation, cell centers (with halo)
#     @inbounds for j ∈ 1:nTh
#         for k ∈ 1:mTh
#             Dhat[k,j] = dudx[k,j+1] - dvdy[k+1,j]
#         end
#     end

#     # very loosely normalizing the inputs to the neural net, not sure if this actually helps/is needed
#     # but I'm leaving the code as it was for all experiments
#     ζ = S.parameters.T.((ζ .- .2) / 15)
#     D = S.parameters.T.(((D .- .3)/ 15))
#     Dhat = S.parameters.T.(((Dhat .- 4e-9) / 10))

#     # move Dhat to cell corners for T12 NN
#     Ixy!(Dhatq, Dhat)

#     # move zeta, D, Dhat to cell centers for T11, T22 NN
#     Ixy!(ζT, ζ)
#     Ixy!(DT, D)
#     Ixy!(DhatT, Dhatq)

#     # Defining two models, the first is temporarily called Su, which will output just T12 and the 
#     # second model, temporarily called Sv, which will output T11 and T22

#     Su_input = Array{T}(undef, nqx, nqy, 3, 1)
#     Su_input[:,:,1,1] .= ζ
#     Su_input[:,:,2,1] .= D
#     Su_input[:,:,3,1] .= Dhatq

#     Sv_input = Array{T}(undef, nx, ny, 3, 1)
#     Sv_input[:,:,1,1] .= ζT
#     Sv_input[:,:,2,1] .= DT
#     Sv_input[:,:,3,1] .= DhatT

#     T12 .= Float64.(Lux.apply(Su_layers, Su_input, model_Su[1], model_Su[2])[1])[:, :, 1, 1]
#     result = Float64.(Lux.apply(Sv_layers, Sv_input, model_Sv[1], model_Sv[2])[1])
#     T11 .= result[:, :, 1, 1]
#     T22 .= result[:, :, 2, 1]

#     ∂x!(dT11dx, T11)
#     ∂y!(dT12dy, T12)

#     ∂x!(dT12dx, T12)
#     ∂y!(dT22dy, T22)

#     @inbounds for j in 1:nuy
#         for k in 1:nux
#             S_u[k,j] = scale * (dT11dx[k,j] + dT12dy[k+1,j])
#         end
#     end

#     @inbounds for j in 1:nvy
#         for k in 1:nvx
#             S_v[k,j] = scale * (dT22dy[k,j] + dT12dx[k,j+1])
#         end
#     end

# end
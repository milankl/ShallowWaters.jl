# juls - A 16bit shallow water model
![sst](figs/sst_posit16.png?raw=true "SST")

A shallow water model written in Julia, which allows arithmetic operations with various number types: 16/32/64bit floats; Arbitrary precision floats (Julia's BigFloat environment); Arbitrary precision posits via the SigmoidNumber package.

Juls is fully-explicit with a Smagorinsky-like biharmonic diffusion operator, the advective terms are written in the vector-invariant form and discretized with either the energy and enstrophy conserving scheme by Arakawa and Hsu, 1990 or the simpler Sadourny 1975 enstrophy conserving scheme. Tracer advection (sofar only passive) is implemented with a semi-Lagrangian advection scheme, that allows for really large time steps, similar to Diamantakis, 2014. Runge Kutte 4th order is used for pressure, advective and coriolis terms and the continuity equation. Semi-implicit time stepping for the diffusive terms (biharmonic diffusion and bottom friction).

Forcing is either with a wind-stress applied to the momentum equations, or surface forcing of the contuinity equation (aka Newtonian Cooling) is possible.Boundary conditions are either periodic (only in x direction) or super-slip/free-slip/partial-slip/no-slip for non-periodic BCs. Output of all prognositc and various diagnostic variables & tendencies is done via NetCDF.

Requires Julia either v0.6, v0.7 or v1.0 and NetCDF.

# HOW TO USE

Change the parameters of your model run in ```parameters.jl``` and then do
```
julia run_juls.jl
```

# THE EQUATIONS

The non-linear shallow water model plus tracer equation is

          ∂u/∂t + (u⃗⋅∇)u - f*v = -g*∂η/∂x - c_D*|u⃗|*u + ∇⋅ν*∇(∇²u) + Fx(x,y)     (1)
          ∂v/∂t + (u⃗⋅∇)v + f*u = -g*∂η/∂y - c_D*|u⃗|*v + ∇⋅ν*∇(∇²v) + Fy(x,y)     (2)
          ∂η/∂t = -∇⋅(u⃗h) + γ*(η_ref - η)                                        (3)
          ∂ϕ/∂t = -u⃗⋅∇ϕ                                                          (4)

with the prognostic variables velocity u⃗ = (u,v) and sea surface heigth η. The layer thickness is h = η + H(x,y). The Coriolis parameter is f = f₀ + βy with beta-plane approximation. The graviational acceleration is g. Bottom friction is either quadratic with drag coefficient c_D or linear with inverse time scale r. Diffusion is realized with a biharmonic diffusion operator, with either a constant viscosity coefficient ν, or a Smagorinsky-like coefficient that scales as ν = c_Smag*|D|, with deformation rate |D| = √((∂u/∂x - ∂v/∂y)² + (∂u/∂y + ∂v/∂x)²). Wind forcing Fx is constant in time, but may vary in space.

The linear shallow water model equivalent is

          ∂u/∂t - f*v = -g*∂η/∂x - r*u + ∇⋅ν*∇(∇²u) + Fx(x,y)     (1)
          ∂v/∂t + f*u = -g*∂η/∂y - r*v + ∇⋅ν*∇(∇²v) + Fy(x,y)     (2)
          ∂η/∂t = -H*∇⋅u⃗ + γ*(η_ref - η)                          (3)
          ∂ϕ/∂t = -u⃗⋅∇ϕ                                           (4)

Juls discretises the equation on an equi-distant Arakawa C-grid, with 2nd order finite-difference operators. Boundary conditions are implemented via a ghost-point copy and each variable has a halo of variable size to account for different stencil sizes of various operators. 

Juls splits the time steps for various terms: Runge Kutta 4th order scheme for the fast varying terms. The diffusive terms (bottom friction and diffusion) are solved semi-implicitly every n-th time step. The tracer equation is solved with a semi-Lagrangian scheme that uses usually much larger time steps.

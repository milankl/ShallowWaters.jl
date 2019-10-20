[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.com/milankl/Juls.jl.svg?branch=master)](https://travis-ci.com/milankl/Juls.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/milankl/Juls.jl?svg=true)](https://ci.appveyor.com/project/milankl/Juls-jl)
[![Build Status](https://api.cirrus-ci.com/github/milankl/Juls.jl.svg)](https://cirrus-ci.com/github/milankl/Juls.jl)


# Juls.jl - A type-flexible 16bit shallow water model
![sst](figs/sst_posit16.png?raw=true "SST")

A shallow water model with a focus on type-flexibility and 16bit number formats. Juls allows for Float64/32/16, BigFloat/[ArbFloat](https://github.com/JeffreySarnoff/ArbNumerics.jl), [Posit32/16](https://github.com/milankl/SoftPosit.jl), [BFloat16](https://github.com/JuliaComputing/BFloat16s.jl), [Sonum16](https://github.com/milankl/Sonums.jl) and in general every number format with arithmetics and conversions implemented.

Juls is fully-explicit with an energy and enstrophy conserving advection scheme and a Smagorinsky-like biharmonic diffusion operator. Tracer advection is implemented with a semi-Lagrangian advection scheme. Runge-Kutta 4th-order is used for pressure, advective and Coriolis terms and the continuity equation. Semi-implicit time stepping for diffusion and bottom friction. Boundary conditions are either periodic (only in x direction) or non-periodic super-slip, free-slip, partial-slip, or no-slip. Output via NetCDF.

# How to use

You find the default parameters in `src/DefaultParameters.jl`. They can be changed with keyword arguments
```julia
julia> Prog = RunJuls(Float32,Ndays=10,g=10,H=500,Fx0=0.12);
Starting Juls on Sun, 20 Oct 2019 19:58:25 without output.
100% Integration done in 4.65s.
```
or by creating a Parameter struct
```julia
julia> P = Parameter(bc="nonperiodic",wind_forcing_x="double_gyre",L_ratio=1,nx=128);
julia> Prog = RunJuls(P);
```

# Installation
```julia
julia> ] add https://github.com/milankl/Juls.jl
```

# THE EQUATIONS

The non-linear shallow water model plus tracer equation is

          ∂u/∂t + (u⃗⋅∇)u - f*v = -g*∂η/∂x - c_D*|u⃗|*u + ∇⋅ν*∇(∇²u) + Fx(x,y)     (1)
          ∂v/∂t + (u⃗⋅∇)v + f*u = -g*∂η/∂y - c_D*|u⃗|*v + ∇⋅ν*∇(∇²v) + Fy(x,y)     (2)
          ∂η/∂t = -∇⋅(u⃗h) + γ*(η_ref - η) + Fηt(t)*Fη(x,y)                       (3)
          ∂ϕ/∂t = -u⃗⋅∇ϕ                                                          (4)

with the prognostic variables velocity u⃗ = (u,v) and sea surface heigth η. The layer thickness is h = η + H(x,y). The Coriolis parameter is f = f₀ + βy with beta-plane approximation. The graviational acceleration is g. Bottom friction is either quadratic with drag coefficient c_D or linear with inverse time scale r. Diffusion is realized with a biharmonic diffusion operator, with either a constant viscosity coefficient ν, or a Smagorinsky-like coefficient that scales as ν = c_Smag*|D|, with deformation rate |D| = √((∂u/∂x - ∂v/∂y)² + (∂u/∂y + ∂v/∂x)²). Wind forcing Fx is constant in time, but may vary in space.

The linear shallow water model equivalent is

          ∂u/∂t - f*v = -g*∂η/∂x - r*u + ∇⋅ν*∇(∇²u) + Fx(x,y)     (1)
          ∂v/∂t + f*u = -g*∂η/∂y - r*v + ∇⋅ν*∇(∇²v) + Fy(x,y)     (2)
          ∂η/∂t = -H*∇⋅u⃗ + γ*(η_ref - η) + Fηt(t)*Fη(x,y)         (3)
          ∂ϕ/∂t = -u⃗⋅∇ϕ                                           (4)

Juls discretises the equation on an equi-distant Arakawa C-grid, with 2nd order finite-difference operators. Boundary conditions are implemented via a ghost-point copy and each variable has a halo of variable size to account for different stencil sizes of various operators.

Juls splits the time steps for various terms: Runge Kutta 4th order scheme for the fast varying terms. The diffusive terms (bottom friction and diffusion) are solved semi-implicitly every n-th time step. The tracer equation is solved with a semi-Lagrangian scheme that uses usually much larger time steps.

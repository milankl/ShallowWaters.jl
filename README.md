# ShallowWaters.jl - A type-flexible 16-bit shallow water model
[![CI](https://github.com/milankl/ShallowWaters.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/milankl/ShallowWaters.jl/actions/workflows/CI.yml)
[![DOI](https://zenodo.org/badge/132787050.svg)](https://zenodo.org/badge/latestdoi/132787050)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://milankl.github.io/ShallowWaters.jl/dev)

![sst](figs/isambard_float16.png?raw=true "Float16 simulation with ShallowWaters.jl on Isambard's A64FX")

A shallow water model with a focus on type-flexibility and 16-bit number formats. ShallowWaters allows for Float64/32/16, 
[Posit32/16/8](https://github.com/milankl/SoftPosit.jl), [BFloat16](https://github.com/JuliaComputing/BFloat16s.jl), 
[LogFixPoint16](https://github.com/milankl/LogFixPoint16s.jl), [Sonum16](https://github.com/milankl/Sonums.jl), 
[Float32/16 & BFloat16 with stochastic rounding](https://github.com/milankl/StochasticRounding.jl) and in 
general every number format with arithmetics and conversions implemented. ShallowWaters also allows for
mixed-precision and reduced precision communication.

ShallowWaters uses an energy and enstrophy conserving advection scheme and a Smagorinsky-like biharmonic diffusion operator. 
Tracer advection is implemented with a semi-Lagrangian advection scheme. Strong stability-preserving Runge-Kutta schemes of
various orders and stages are used with a semi-implicit treatment of the continuity equation. Boundary conditions are either 
periodic (only in x direction) or non-periodic super-slip, free-slip, partial-slip, or no-slip.
Output via [NetCDF](https://github.com/JuliaGeo/NetCDF.jl).

Please feel free to raise an [issue](https://github.com/milankl/ShallowWaters.jl/issues) if you discover bugs or have an idea how to improve ShallowWaters.

Requires: Julia 1.2 or higher

## How to use

```julia
julia> using ShallowWaters
julia> run_model()
Starting ShallowWaters on Thu, 20 Jan 2022 15:44:30 without output.
100% Integration done in 0.61s.
```
You just successfully ran ShallowWaters.jl! For more examples and arguments to pass on to `run_model` see the [documentation](https://milankl.github.io/ShallowWaters.jl/dev).

## Installation

ShallowWaters.jl is a registered package, so simply do

```julia
julia> ] add ShallowWaters
```

## References

ShallowWaters.jl was used and is described in more detail in  

>[1] Klöwer M, PD Düben, and TN Palmer, 2020. Number formats, error mitigation and scope for 16-bit arithmetics in weather and climate modelling analysed with a shallow water model. __Journal of Advances in Modeling Earth Systems__, [10.1029/2020MS002246](https://dx.doi.org/10.1029/2020MS002246)

>[2] Klöwer M, PD Düben, and TN Palmer, 2019. Posits as an alternative to floats for weather and climate models. In: __Proceedings of the Conference for Next Generation Arithmetic 2019__, Singapore, __ACM__, [10.1145/3316279.3316281](https://dx.doi.org/10.1145/3316279.3316281)

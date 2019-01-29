# juls - A 16bit shallow water model
![sst](figs/sst_posit16.png?raw=true "SST")

A shallow water model (but hopefully going to be 2-Layer primitive equation model) written in Julia, with a focus on flexible number types: 16/32/64bit floats; Arbitrary precision floats (Julia's BigFloat environment); Arbitrary precision posits via the SigmoidNumber package; Maybe in the future also Integers.

Juls is fully-explicit with a Smagorinsky-like biharmonic diffusion operator, the advective terms are written in the vector-invariant form and discretized with either the energy and enstrophy conserving scheme by Arakawa and Hsu, 1990 or the simpler Sadourny 1975 enstrophy conserving scheme. Tracer advection (sofar only passive) is implemented with a semi-Lagrangian advection scheme, that allows for really large time steps, similar to Diamantakis, 2014. Runge Kutte 4th order is used for pressure, advective and coriolis terms and the continuity equation. Semi-implicit time stepping for the diffusive terms (biharmonic diffusion and bottom friction).

Forcing is either with a wind-stress applied to the momentum equations, or surface forcing of the contuinity equation (aka Newtonian Cooling) is possible.

Boundary conditions are either periodic (only in x direction) or super-slip/free-slip/partial-slip/no-slip for non-periodic BCs.

Output is done via NetCDF.

Requires Julia either v0.6, v0.7 or v1.0 and NetCDF.

# TODO

- from one layer to two (or n-) layers
- split external mode (gravity waves) from internal modes - a potential speed up of 10-100x
- parallelize via domain decomposition (all operators are written in a generic way, applicable to every subdomain too, such that only the ghost point copy needs to include communication between threads)
- split the tracer relaxation into a tracer source and a tracer sink with different time scales
- include an atmospheric setting
- bicubic interpolation for tracer advection

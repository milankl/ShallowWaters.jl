# juls
A shallow water model (but going to be 2-Layer primitive equation model) written in Julia, with a focus on flexible number types: 16/32/64bit floats; Arbitrary precision floats (Julia's BigFloat environment); Arbitrary precision posits via the SigmoidNumber package (maybe in the future also Integers).

Juls is fully-explicit with a Smagorinsky-like biharmonic diffusion operator, the advective terms are written in the vector-invariant form and discretized with either the energy and enstrophy conserving scheme by Arakawa and Hsu, 1990 or the simpler Sadourny 1975 enstrophy conserving scheme. Tracer advection (sofar only passive) is implemented with a semi-Lagrangian advection scheme, that allows for really large time steps, similar to Diamantakis, 2014.

Forcing is either with a wind-stress applied to the momentum equations, or surface forcing of the contuinity equation (aka Newtonian Cooling) is possible.

Boundary conditions are either periodic (only in x direction) or super-slip/free-slip/partial-slip/no-slip for non-periodic BCs.

Output is done via NetCDF.

Requires Julia v0.6 or v0.7 and NetCDF.

# TODO

- from one layer to two (or n-) layers
- split external mode (gravity waves) from internal modes - a potential speed up of 100x
- parallelize via domain decomposition (all operators are written in a generic way, applicable to every subdomain too, such that only the ghost point copy needs to include communication between threads)
- split the tracer relaxation into a tracer source and a tracer sink with different time scales

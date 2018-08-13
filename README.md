# juls
A (going to be) 2-Layer primitive equation model written in Julia, with a focus on flexible number types: 16/32/64bit floats; Arbitrary precision floats (Julia's BigFloat environment); Arbitrary precision posits via the SigmoidNumber package; maybe in the future also Integers.

Juls is fully-explicit with a Smagorinsky-like biharmonic diffusion operator, the advective terms are written in the vector-invariant form and discretized with the Sadourny 1975 enstrophy conserving scheme (to be changed to Arakawa and Hsu, 1990 scheme).

Boundary conditions are either periodic (only in x direction) or super-slip/free-slip/partial-slip/no-slip for non-periodic BCs.

Output is done via NetCDF.

Requires Julia v0.6, NetCDF.

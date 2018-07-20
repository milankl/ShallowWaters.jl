# juls
A (sofar only) 1-Layer shallow water model written in Julia, with a focus on flexible number types: 64bit/32bit Floats; Arbitrary precision floats (Julia's BigFloat environment); Arbitrary preicision Posits via the SigmoidNumber package; Integers.

Juls is fully-explicit with a Smagorinsky-like biharmonic diffusion operator, the advective terms are written in the vector-invariant form and discretized with the Sadourny 1975 enstrophy conserving scheme (to be changed to Arakawa and Hsu, 1990 scheme).

Boundary conditions are either periodic (only in x direction) or super-slip/free-slip/partial-slip/no-slip for non-periodic BCs.

Output is done via NetCDF.

Requires Julia v0.6, NetCDF.

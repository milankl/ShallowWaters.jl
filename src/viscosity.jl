# the dx^2 factor is removed from the Laplace operators
# ν can be used for harmonic as well as for biharmonic
#const ν = Numtype(540*128/nx/Δ^2)

# Smagorinsky-like biharmonic diffusion
const cSmag = Numtype(-c_smag)
const νA = 1000.0/Δ^2

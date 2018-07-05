# the dx^2 factor is removed from the Laplace operators
# ν can be used for harmonic as well as for biharmonic
const ν = Numtype(540*128/nx/dx^2)

# Smagorinsky-like biharmonic diffusion
const cSmag = Numtype(0.01)

#=
LLu_Smag = dx[ νSmag dx(L(u))] + dy[ νSmag dy(L(u))]
LLv_Smag = dx[ νSmag dx(L(v))] + dy[ νSmag dy(L(v))]
=#

#= Smagorinsky-like biharmonic viscosity

    Viscosity = ∇ ⋅ (cSmag Δ⁴ |D| ∇∇² ⃗u)

The Δ⁴-scaling is omitted as gradient operators are dimensionless. =#

const cSmag = Numtype(-c_smag)

#= Constant viscosity coefficient that can be used for harmonic or biharmonic
diffusion operators

    Viscosity = ν∇² ⃗u      or  = ν∇⁴ ⃗u
=#

# const ν = Numtype(540*128/nx/Δ^2)

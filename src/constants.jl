# Set up constants to be used within the RHS rhs.jl
# Set typeof according to the used number type Numtype

# for the ghost point copy/tangential boundary conditions
const one_minus_α = Numtype(1-lbc)

# for the laplace operator
const minus_4 = Numtype(-4.)

# for the interpolation functions
const one_half = Numtype(0.5)
const one_quarter = Numtype(0.25)

# will be used for the Bernoulli potential
const g = Numtype(gravity)

# currently not used
# const zeero = Numtype(0.)
# const α = Numtype(lbc)
# const minus_α = Numtype(-lbc)
# const one_minus_α_half = Numtype(1-0.5*lbc)
# const α_over_Δ = Numtype(lbc/Δ)
# const minus_3_minus_α = Numtype(-3-lbc)

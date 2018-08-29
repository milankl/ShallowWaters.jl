# Set up constants to be used within the RHS rhs.jl
# Set typeof according to the used number type Numtype

# for the ghost point copy/tangential boundary conditions
const one_minus_α = Numtype(1-lbc)

# for the laplace operator
const minus_4 = Numtype(-4.)

# for the interpolation functions
const one_half = Numtype(0.5)
const one_twelve = Numtype(1/12)
const one_quarter = Numtype(0.25)

# will be used for the Bernoulli potential
const g = Numtype(gravity)

# for the bottom friction
const c_D = Numtype(drag)

# frequency [1/s] of interface relaxation (for non-dimensional gradients, γ contains the grid spacing Δ)
const γ = Numtype(Δ/(t_relax*3600*24))

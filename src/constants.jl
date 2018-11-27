# Set up constants to be used within the RHS rhs.jl
# Set typeof according to the used number type Numtype

# Runge-Kutta 3rd/4th order coefficients
if RKo == 3     # version 2
    const RKa = Numtype.([1/4,0.,3/4])
    const RKb = Numtype.([1/3,2/3])
elseif RKo == 4
    const RKa = Numtype.([1/6,1/3,1/3,1/6])
    const RKb = Numtype.([.5,.5,1.])
end

# for the ghost point copy/tangential boundary conditions
const one_minus_α = Numtype(1-lbc)

# for the laplace operator
const minus_4 = Numtype(-4.)

# for the interpolation functions
const one_half = Numtype(0.5)
const one_twelve = Numtype(1/12)
const one_quart = Numtype(0.25)

# will be used for the Bernoulli potential
const g = Numtype(gravity)

# for the bottom friction (include grid spacing as gradient operators are dimensionless)
const c_D = Numtype(Δ*drag)
const r_B = Numtype(Δ/(τdrag*24*3600))  # [1/s]

# frequency [1/s] of interface relaxation (for non-dimensional gradients, γ contains the grid spacing Δ)
const γ = Numtype(Δ/(t_relax*3600*24))

# for biharmonic diffusion
const cSmag = Numtype(-c_smag)
const νB = Numtype(-ν_const/30000)   # linear scaling based on 540m^s/s at Δ=30km

# for semi-Lagrangian advection / interpolation
const zeero = Numtype(0.)
const oone = Numtype(1.)
const r_SST = Numtype(dtadvint/(τSST*3600*24))    # [1]

# for analysis the old operators with boundary conditions are kept. Constants for these
const one_over_Δ = Numtype(1/Δ)
const α = Numtype(lbc)
const minus_α = Numtype(-lbc)
const one_minus_α_half = Numtype(1-0.5*lbc)
const α_over_Δ = Numtype(lbc/Δ)
const one_quarter = Numtype(0.25)
const minus_3_minus_α = Numtype(-3-lbc)

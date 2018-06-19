# define as constants what is needed inside the rhs

const nx = 64
const Lx = 512e3
const L_ratio = 3

const g = 10.
const H = 500.

const cfl = 0.9
const Ndays = 10

# boundary condtions
const bc_x = "periodic"    # or "free-slip" or "no-slip"
const bc_y = "free-slip"   # or "no-slip"

const run_id = 0           #TODO make automatic
const c_D = 1e-5

const output = 1

const lat_0 = 30.

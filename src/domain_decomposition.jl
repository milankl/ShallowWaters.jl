""" Factorizes the positive Integer n into its two largest factors.
e.g. factorization(12) yields (3,4) as 3*4=12 (instead of 2*6)."""
function factorization(n::Int)
    s = floor(Int, sqrt(n)))
    m = mod(n,s)

    while m != 0
        s -= 1
        m = mod(n,s)
    end

    # Put more processes in x or y direction which has more grid cells.
    # This attempts to make the subdomains as quadratic as possible.
    if L_ratio <= 1
        return s,n ÷ s
    else
        return n ÷ s,s
    end
end

""" Returns a matrix of a cartesian grid domain decomposition
to determine the ranks of the neighbouring processes.
Ranks start from 0.

E.g. neighours_mat(12) returns a 3x4 matrix
0  3  6   9
1  4  7  10
2  5  8  11

as (3,4) are the two largest factors of 12."""
function neighbours_mat(nprocs::Int)
    f1,f2 = factorization(nprocs)
    return Array(reshape(0:nprocs-1,f1,f2))
end

""" Returns an array of [n,e,s,w] (North, East, South, West) neighbours
for a given rank in domain decomposition."""
function neighbours(prank::Int,nprocs::Int)
    M = neighbours_mat(nprocs)
    m,n = size(M)

    prank >= 0 && prank < nprocs || throw(error("Processor rank must be within (0,$(nprocs-1)), $prank was given."))

    # get indices of process rank
    i,j = mod(prank,m)+1,(prank ÷ m)+1

    # find neighbours, -1 for non existent
    # nnb north neighbour, snb south neighbour etc
    nnb = if (j+1 > n); -1 else M[i,j+1] end
    snb = if (j-1 < 1); -1 else M[i,j-1] end
    if bc_x == "periodic"
        enb = M[mod(i,m)+1,j]
        wnb = M[mod(i-2,m)+1,j]
    else
        enb = if (i+1 > m); -1 else M[i+1,j] end
        wnb = if (i-1 < 1); -1 else M[i-1,j] end
    end
    return [nnb,enb,snb,wnb],M
end

function subdomain_grid(M::AbstractMatrix,nbours::AbstractVector,prank::Int,nprocs::Int)

    m,n = size(M)

    # throw error if domain is not equally distributable among processes
    mod(nx,m) == 0 || mod(ny,n) == 0 ||
        throw(DimensionMismatch("Domainsize $nx,$ny is incompatible with domain decomposition ($m,$n)."))

    sub_dom = Dict()

    # number of grid points on the T-grid for given direction x,y
    sub_dom["nx"] = nx ÷ m
    sub_dom["ny"] = ny ÷ n

    # number of grid points on the u-grid for given direction x,y
    if nbours[4] != -1  # if west neighbour exists
        sub_dom["nux"] = sub_dom["nx"]
    else
        sub_dom["nux"] = sub_dom["nx"]-1
    end

    sub_dom["nuy"] = sub_dom["ny"]

    # number of grid points on the v-grid for given direction x,y
    sub_dom["nvx"] = sub_dom["nx"]

    if nbours[1] == -1 # if north neighbour does not exist
        sub_dom["nvy"] = sub_dom["ny"]-1
    else
        sub_dom["nvy"] = sub_dom["ny"]
    end

    # number of grid points on the q-grid for given direction x,y
    if nbours[2] == -1 # if east neighbour does not exit
        sub_dom["nqx"] = sub_dom["nx"]+1
    else
        sub_dom["nqx"] = sub_dom["nx"]
    end

    if nbours[3] == -1 # if south neighbour does not exist
        sub_dom["nqy"] = sub_dom["ny"]+1
    else
        sub_dom["nqy"] = sub_dom["ny"]
    end

    # total number of T,u,v,q-points
    sub_dom["nT"] = sub_dom["nx"]*sub_dom["ny"]
    sub_dom["nu"] = sub_dom["nux"]*sub_dom["nuy"]
    sub_dom["nv"] = sub_dom["nvx"]*sub_dom["nvy"]
    sub_dom["nq"] = sub_dom["nqx"]*sub_dom["nqy"]

    return sub_dom
end

# MPI.Init()
# const comm = MPI.COMM_WORLD
#
# const proc_rank = MPI.Comm_rank(comm)
# const nprocs = MPI.Comm_size(comm)
#
# const nbours,M = neighbours(proc_rank,nprocs)
# const sub_dom = subdomain_grid(M,nbours,proc_rank,nprocs)

# for testing
#nbours,M = neighbours(0,1)
#sub_dom = subdomain_grid(M,nbours,0,1)

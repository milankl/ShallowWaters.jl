""" Factorizes the positive Integer n into its two largest factors.
e.g. factorization(12) yields (3,4) as 3*4=12 (instead of 2*6)."""
function factorization(n::Int)
    s = Int(floor(sqrt(n)))
    m = mod(n,s)

    while m != 0
        s -= 1
        m = mod(n,s)
    end

    return s,n ÷ s
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
    return Array(reshape(0:size-1,f1,f2))
end

""" Returns an array of [n,e,s,w] (North, East, South, West) neighbours
for a given rank in domain decomposition."""
function neighbours(prank::Int,nprocs::Int)
    M = neighbours_mat(nprocs)
    m,n = size(M)

    # get indices of process rank
    i,j = mod(prank,m)+1,(prank ÷ m)+1

    # find neighbours, -1 for non existent
    # nnb north neighbour, snb south neighbour etc
    nnb = if (j+1 > n); -1 else M[i,j+1] end
    snb = if (j-1 < 1); -1 else M[i,j-1] end
    if bc_x == "periodic"
        enb = M[mod(i+1,m),j]
        wnb = M[mod(i-1,m),j]
    else
        enb = if (i+1 > m); -1 else M[i+1,j] end
        wnb = if (i-1 < 1); -1 else M[i-1,j] end
    end
    return [nnb,enb,snb,wnb]
end

function AHα!(α::AbstractMatrix,q::AbstractMatrix)
    #= Linear combination α of potential voriticity q according
    to the energy and enstrophy conserving scheme of Arakawa and Hsu, 1990 =#

    m,n = size(α)
    @boundscheck (m+1,n+1) == size(q) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            α[i,j] = one_twelve*(q[i,j] + q[i,j+1] + q[i+1,j+1])
        end
    end
end

function AHβ!(β::AbstractMatrix,q::AbstractMatrix)
    #= Linear combination δ of potential voriticity q according
    to the energy and enstrophy conserving scheme of Arakawa and Hsu, 1990 =#

    m,n = size(β)
    @boundscheck (m,n+1) == size(q) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        β[1,j] = one_twelve*(q[1,j] + q[1,j+1] + q[end,j+1])
        for i ∈ 2:m
            β[i,j] = one_twelve*(q[i,j] + q[i,j+1] + q[i-1,j+1])
        end
    end
end

function AHγ!(γ::AbstractMatrix,q::AbstractMatrix)
    #= Linear combination γ of potential voriticity q according
    to the energy and enstrophy conserving scheme of Arakawa and Hsu, 1990 =#

    m,n = size(γ)
    @boundscheck (m,n+1) == size(q) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        γ[1,j] = one_twelve*(q[1,j] + q[1,j+1] + q[end,j])
        for i ∈ 2:m
            γ[i,j] = one_twelve*(q[i,j] + q[i-1,j] + q[i,j+1])
        end
    end
end

function AHδ!(δ::AbstractMatrix,q::AbstractMatrix)
    #= Linear combination β of potential voriticity q according
    to the energy and enstrophy conserving scheme of Arakawa and Hsu, 1990 =#

    m,n = size(δ)
    @boundscheck (m+1,n+1) == size(q) || throw(BoundsError())

    @inbounds for j ∈ 1:n
        for i ∈ 1:m
            δ[i,j] = one_twelve*(q[i,j] + q[i,j+1] + q[i+1,j])
        end
    end
end

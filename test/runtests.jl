using Juls
using Test

@testset "StartFromRest"
    for T in (Float32,Float64)
        u,v,η = RunJuls(T,output=false)
        @test all(u .== zero(T))
        @test all(v .== zero(T))
        @test all(η .== zero(T))
    end
end

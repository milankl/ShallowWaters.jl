using Juls
using Test

@testset "NoForcing" begin
    Prog = RunJuls(Ndays=1,Fx0=0)
    @test all(Prog.u .== 0.0f0)
    @test all(Prog.v .== 0.0f0)
    @test all(Prog.η .== 0.0f0)
end

@testset "BoundaryConditions" begin
    Prog = RunJuls(Ndays=1,bc="periodic")
    @test all(abs(Prog.η) .< 10)    # sea surface height shouldn't exceed +-10m

    Prog = RunJuls(Ndays=1,bc="nonperiodic")
    @test all(abs(Prog.η) .< 10)

    Prog = RunJuls(Ndays=1,α=0)     # free-slip
    @test all(abs(Prog.η) .< 10)
end

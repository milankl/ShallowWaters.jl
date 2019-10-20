using Juls
using Test

@testset "NoForcing" begin
    Prog = RunJuls(Ndays=1,Fx0=0)
    @test all(Prog.u .== 0.0f0)
    @test all(Prog.v .== 0.0f0)
    @test all(Prog.η .== 0.0f0)
end

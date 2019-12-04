#using Juls
using Test

@testset "No Forcing" begin
    Prog = RunJuls(Ndays=1,Fx0=0)
    @test all(Prog.u .== 0.0f0)
    @test all(Prog.v .== 0.0f0)
    @test all(Prog.η .== 0.0f0)
end

@testset "Boundary Conditions" begin
    Prog = RunJuls(Ndays=1,bc="periodic")
    @test all(abs.(Prog.η) .< 10)    # sea surface height shouldn't exceed +-10m

    Prog = RunJuls(Ndays=1,bc="nonperiodic")
    @test all(abs.(Prog.η) .< 10)

    Prog = RunJuls(Ndays=1,α=0)     # free-slip
    @test all(abs.(Prog.η) .< 10)
end

@testset "Continuity Forcing" begin
    Prog = RunJuls(Ndays=1,surface_relax=true)
    @test all(abs.(Prog.η) .< 10)    # sea surface height shouldn't exceed +-10m
    @test all(Prog.u .!= 0.0f0)

    Prog = RunJuls(Ndays=1,surface_forcing=true)
    @test all(abs.(Prog.η) .< 10)
    @test all(Prog.u .!= 0.0f0)
end

@testset "Mixed Precision" begin
    Prog = RunJuls(Float32,Tprog=Float64,Ndays=1)
    @test all(abs.(Prog.η) .< 10)    # sea surface height shouldn't exceed +-10m
    @test all(Prog.u .!= 0.0f0)

    Prog = RunJuls(Float16,Tprog=Float32,Ndays=1)
    @test all(abs.(Prog.η) .< 10)    # sea surface height shouldn't exceed +-10m
    @test all(Prog.u .!= 0.0f0)
end

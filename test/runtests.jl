using ShallowWaters
using Test

@testset "No Forcing" begin
    Prog = run_model(Fx0=0)
    @test all(Prog.u .== 0.0f0)
    @test all(Prog.v .== 0.0f0)
    @test all(Prog.η .== 0.0f0)
end

@testset "Boundary Conditions" begin
    Prog = run_model(bc="periodic")
    @test all(abs.(Prog.u) .< 10)    # velocity shouldn't exceed +-10m/s

    Prog = run_model(bc="nonperiodic")
    @test all(abs.(Prog.u) .< 10)

    Prog = run_model(α=0)     # free-slip
    @test all(abs.(Prog.u) .< 10)
end

@testset "Linear dynamics" begin
    P = run_model(dynamics="linear")
    @test all(abs.(P.u) .< 10)   # velocity shouldn't exceed +-10m/s
end

@testset "Continuity Forcing" begin
    Prog = run_model(surface_relax=true)
    @test all(abs.(Prog.u) .< 10)    # velocity shouldn't exceed +-10m/s
    @test all(Prog.u .!= 0.0f0)

    Prog = run_model(surface_forcing=true)
    @test all(abs.(Prog.u) .< 10)
    @test all(Prog.u .!= 0.0f0)
end

@testset "Mixed Precision" begin
    Prog = run_model(Float32,Tprog=Float64)
    @test all(abs.(Prog.u) .< 10)    # velocity shouldn't exceed +-10m/s
    @test all(Prog.u .!= 0.0f0)

    Prog = run_model(Float16,Tprog=Float32,Ndays=50)
    @test all(abs.(Prog.u) .< 10)    # velocity shouldn't exceed +-10m/s
    @test all(Prog.u .!= 0.0f0)
end

@testset "RK3 tests" begin
    Prog = run_model(time_scheme="RK",RKo=3,cfl=0.6)
    @test all(abs.(Prog.u) .< 10)    # velocity shouldn't exceed +-10m/s
    @test all(Prog.u .!= 0.0f0)
end

@testset "SSPRK2 tests" begin
    Ndays = 500
    Prog = run_model(time_scheme="SSPRK2",RKs=2,cfl=0.7;Ndays)
    @test all(abs.(Prog.u) .< 10)    # velocity shouldn't exceed +-10m/s
    @test all(Prog.u .!= 0.0f0)

    Prog = run_model(time_scheme="SSPRK2",RKs=3,cfl=1.2;Ndays)
    @test all(abs.(Prog.u) .< 10)    # velocity shouldn't exceed +-10m/s
    @test all(Prog.u .!= 0.0f0)

    Prog = run_model(time_scheme="SSPRK2",RKs=4,cfl=1.3;Ndays)
    @test all(abs.(Prog.u) .< 10)    # velocity shouldn't exceed +-10m/s
    @test all(Prog.u .!= 0.0f0)
end

@testset "SSPRK3 tests" begin
    Ndays = 500
    Prog = run_model(time_scheme="SSPRK3",RKn=2,cfl=1.2;Ndays)
    @test all(abs.(Prog.u) .< 10)    # velocity shouldn't exceed +-10m/s
    @test all(Prog.u .!= 0.0f0)

    Prog = run_model(time_scheme="SSPRK3",RKn=3,cfl=2.5;Ndays)
    @test all(abs.(Prog.u) .< 10)    # velocity shouldn't exceed +-10m/s
    @test all(Prog.u .!= 0.0f0)
end

@testset "Scaling" begin
    P64 = run_model(nx=100,scale=64);
    P1 = run_model(nx=100,scale=1);

    @test P64.u == P1.u
    @test P64.v == P1.v
    @test P64.η == P1.η
    @test P64.sst == P1.sst
end
using AnalyticMHDTestSolutions, Test

@testset "Riemann Sod-Shock" begin

    par = RiemannParameters(Ul=100.0, Mach=10.0, t=1.5)

    sol = solve([86.0], par)

    @test sol.rho3 ≈ 0.37667409437005994
    @test sol.rho4 ≈ 0.4854368932038835
    @test sol.P34  ≈ 13.097397114653086

end

@testset "Riemann CR-Shock" begin

    par = RiemannParameters(Pl=63.400, Pr=0.1, Mach=10.0, dsa_model=4, t=1.5)

    sol = solve([86.0], par)

    @test sol.rho3    ≈ 0.3699882303652922
    @test sol.rho4    ≈ 0.5925766595991485
    @test sol.P34_tot ≈ 12.089335761741005
    @test sol.P4_cr   ≈ 3.583255961157783

end

@testset "DSA models" begin

    #@test_warn "Invalid DSA model selection!" RiemannParameters( Ul = 100.0, Ur = 0.1, dsa_model=5, t = 1.5)

    # KR07
    par = RiemannParameters( Ul = 100.0, Ur = 0.1, dsa_model=0, t = 1.5)
    @test par.acc_function(5.0) ≈ 0.25185919999999995
    # KR13
    par = RiemannParameters( Ul = 100.0, Ur = 0.1, dsa_model=1, t = 1.5)
    @test par.acc_function( 5.0) ≈ 0.09999999999999998
    @test par.acc_function(10.0) ≈ 0.19631644350722818
    @test par.acc_function(25.0) ≈ 0.21152
    # Ryu+19
    par = RiemannParameters( Ul = 100.0, Ur = 0.1, dsa_model=2, t = 1.5)
    @test par.acc_function( 5.0) ≈ 0.017286554080677037
    @test par.acc_function(55.0) ≈ 0.0348
    # CS14
    par = RiemannParameters( Ul = 100.0, Ur = 0.1, dsa_model=3, t = 1.5)
    @test par.acc_function(5.0) ≈ 0.04999999999999999
    # Pfrommer+16
    par = RiemannParameters( Ul = 100.0, Ur = 0.1, dsa_model=4, t = 1.5)
    @test par.acc_function(5.0) == 0.5
end
using AnalyticMHDTestSolutions, Test

@testset "Shocktubes" begin
    
    @testset "Riemann Sod-Shock" begin

        par = RiemannParameters(Ul=100.0, Mach=10.0, t=1.5)

        x = collect(50.0:0.01:100.0)

        sol = solve(x, par)

        @test sol.rho3 ≈ 0.37667409437005994
        @test sol.rho4 ≈ 0.4854368932038835
        @test sol.P34  ≈ 13.097397114653086

        # parameter tests
        @test_nowarn par = RiemannParameters(Pl=100.0, Mach=10.0, t=1.5)
        #@test_nowarn par = RiemannParameters(Pr=1.0, Mach=10.0, t=1.5)
        #@test_nowarn par = RiemannParameters(Ur=1.0, Mach=10.0, t=1.5)

        #@test_throws ErrorException("No initial Pressure or energy values given!") par = RiemannParameters(Mach=10.0, t=1.5)

        #@test_nowarn par = RiemannParameters(Pr=1.0, Mach=10.0, t=1.5)
    
        #@test_throws ErrorException("Ur, Pr and Mach are zero! Can't find solution!") par = RiemannParameters(Ul=100.0, t=1.5)
        #@test_throws ErrorException("Ul, Pl and Mach are zero! Can't find solution!") par = RiemannParameters(Ur=1.0, t=1.5)

        @test_nowarn par = RiemannParameters(Ul=100.0, Ur=1.0, t=1.5)

    end

    @testset "Riemann CR-Shock" begin

        par = RiemannParameters(Pl=63.400, Pr=0.1, Mach=10.0, dsa_model=4, t=1.5)

        x = collect(50.0:0.01:100.0)

        sol = solve(x, par)

        @test sol.rho3    ≈ 0.3699882303652922
        @test sol.rho4    ≈ 0.5925766595991485
        @test sol.P34_tot ≈ 12.089335761741005
        @test sol.P4_cr   ≈ 3.583255961157783

    end

end

@testset "Sedov" begin
    
    @testset "SedovParameters" begin
        
        sedov_par = SedovParameters(0.05, Ndim=2)

        @test sedov_par.rho_s ≈ 4.0
        @test sedov_par.cs_out ≈ 0.10540925533894599

    end

    @testset "Solve Sedov" begin
        
        sedov_par = SedovParameters(0.05, Ndim=2)

        r = collect(0:1.e-6:1)
        sedov_solution = solve(r, sedov_par)

        @test sedov_solution.r_shock ≈ 0.2579946323580024
        @test sedov_solution.v_shock ≈ 1.9349597426850176

        @test maximum(sedov_solution.rho) ≈ 3.999919116451743
        @test maximum(sedov_solution.P) ≈ 4.992021919033185
        @test maximum(sedov_solution.vr) ≈ 1.9349514429992791

    end

end
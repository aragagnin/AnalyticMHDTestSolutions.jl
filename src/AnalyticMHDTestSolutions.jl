module AnalyticMHDTestSolutions

    using DSAModels 
    
    # riemann solvers
    include("sod_shocks/cr_sod_shock_main.jl")
    include("sod_shocks/setup_riemann_parameters.jl")

    # sedov solution
    include("sedov/sedov_structs.jl")
    include("sedov/sedov_fit.jl")
    include("sedov/sedov_solution.jl")

    export RiemannParameters,    # helper function to set up solution
           solve,                # overloaded function to solve riemann problems
           find_xs_first_guess,  # helper function to find initial guess for shock compression

           SedovData,
           SedovParameters,
           SedovSolution,
           get_sedov_solution,    # function to get sedov data and ideal solution from sedovdata
           get_sedov_data_from_gadget,
           get_sedov_solution_from_gadget


    """
    Multiple dispatch for solve function
    """
    # Pure hydro Sod-shock
    """
        solve(x::Array{Float,1}, par::SodParameters)

    Solves a standard Sod shock problem.
    """
    solve(x::Array{Float64,1}, par::SodParameters) = solveSodShock(x, par=par)
    solve(x::Array{Float32,1}, par::SodParameters) = solveSodShock(Float64.(x), par=par)

    # CR Sod shock
    """
        solve(x::Array{Float,1}, par::SodCRParameters_noCRs)

    Solves a Sod shock with cosmic ray acceleration, but without a pre-existing CR component.
    """
    solve(x::Array{Float64,1}, par::SodCRParameters_noCRs)   = solveSodShockCR_noPrepopulation(x, par=par)
    solve(x::Array{Float32,1}, par::SodCRParameters_noCRs)   = solveSodShockCR_noPrepopulation(Float64.(x), par=par)

    """
        solve(x::Array{Float,1}, par::SodCRParameters_withCRs)

    Solves a Sod shock with cosmic ray acceleration without a pre-existing CR component.
    """
    solve(x::Array{Float64,1}, par::SodCRParameters_withCRs) = solveSodShockCR_withPrepopulation(x, par=par)
    solve(x::Array{Float32,1}, par::SodCRParameters_withCRs) = solveSodShockCR_withPrepopulation(Float64.(x), par=par)



end # module

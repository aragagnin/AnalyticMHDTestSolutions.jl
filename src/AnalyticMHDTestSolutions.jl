module AnalyticMHDTestSolutions

    # riemann solvers
    include("cr_dsa_models.jl")
    include("cr_sod_shock_main.jl")
    include("setup_riemann_parameters.jl")
    #include(joinpath(dirname(@__FILE__), "ideal_solutions", "cr_sod_shock.jl"))

    # sedov solution
    include("sedov_structs.jl")
    include("sedov_fit.jl")
    include("sedov_solution.jl")

    export RiemannParameters,    # helper function to set up solution
           solve,                # overloaded function to solve riemann problems
           find_xs_first_guess,  # helper function to find initial guess for shock compression

           SedovData,
           SedovParameters,
           SedovSolution,
           get_sedov_solution,    # function to get sedov data and ideal solution from sedovdata
           get_sedov_data_from_gadget,
           get_sedov_solution_from_gadget



end # module

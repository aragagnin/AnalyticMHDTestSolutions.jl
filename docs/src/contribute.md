```@meta
CurrentModule = AnalyticMHDTestSolutions
DocTestSetup = quote
    using AnalyticMHDTestSolutions
end
```

# Contributing

If you want to contribute to this you are very welcome to do so!

I only ask you to follow some basic procedures:

1. Fork the repository
2. Issue pull-requests from your fork to the main repository
3. Write unit-tests for your newly added solvers
4. Write documentation for your newly added solvers

As a general style idea I would suggest that you write a function `solve` that uses multiple dispatch to handle your newly implemented solver and returns a `struct` with the solution. Something like:

```julia
solve( pos::Array{<:Real}, par::YourSolverParameters ) = your_solver_function( pos, par )
```

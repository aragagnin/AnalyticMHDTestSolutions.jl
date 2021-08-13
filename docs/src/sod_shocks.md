```@meta
CurrentModule = AnalyticMHDTestSolutions
DocTestSetup = quote
    using AnalyticMHDTestSolutions
end
```

# Sod Shocks



## Hydrodynamic Sod Shocks

To get the exact solution to a Sod shock you first need to set up the inital conditions.
You can do this with the helper function [`RiemannParameters`](@ref) that contains all parameters for all possible configurations:

```julia
RiemannParameters( ;rhol::Float64=1.0, rhor::Float64=0.125,      # density left and right (L&R)
                    Pl::Float64=0.0,   Pr::Float64=0.0,          # pressure L&R
                    Ul::Float64=0.0,   Ur::Float64=0.0,          # internal energy L&R
                    P_cr_l::Float64=0.0, P_cr_r::Float64=0.0,    # CR pressure L&R
                    E_cr_l::Float64=0.0, E_cr_r::Float64=0.0,    # CR energy L&R
                    Bl::Array{Float64,1} = zeros(3),             # B-field left
                    Br::Array{Float64,1} = zeros(3),             # B-field right
                    Mach::Float64=0.0,                           # target Mach number
                    t::Float64,                                  # time of the solution
                    x_contact::Float64=70.0,                     # position of the contact discontinuity along the tube
                    γ_th::Float64=5.0/3.0,                       # adiabatic index of the gas
                    γ_cr::Float64=4.0/3.0,                       # adiabatic index of CRs
                    Pe_ratio::Float64=0.01,                      # ratio of proton to electron energy in acceleration
                    thetaB::Float64=0.0,                         # angle between magnetic field and shock normal
                    theta_crit::Float64=(π/4.0),                 # critical angle for B/Shock angle efficiency
                    dsa_model::Int64=-1,                         # diffuse shock acceleration model
                    xs_first_guess::Float64=4.7,                 # first guess of the resulting shock compression
                    verbose::Bool=false)
```

To set up a standard Sod shock you need to supply it with pressure/energy values for left and right state, or with pressure/energy values for the left state and a target Mach number.

A minimal working version would be, for a shock with Mach 10, at time = 1.5:

```julia
par = RiemannParameters(Ul=100.0, Mach=10.0, t=1.5)
```

This returns a parameter object for a pure hydro Sod shock:

```julia
mutable struct SodParameters

    rhol::Float64           # denisty left
    rhor::Float64           # density right
    Pl::Float64             # pressure left
    Pr::Float64             # pressure right
    Ul::Float64             # internal energy left
    Ur::Float64             # internal energy right
    cl::Float64             # soundspeed left
    cr::Float64             # soundspeed right
    M::Float64              # Mach number
    t::Float64              # time
    x_contact::Float64      # position of the contact discontinuity along the tube
    γ_th::Float64           # adiabatic index of the gas
    γ_exp::Float64          # helper variable
    η2::Float64             # helper variable

end
```

## Sod Shock with CR acceleration

You can solve a sod shock with CR acceleration by chosing the desired acceleration function with the keyword argument `dsa_model`.
Possible choices are:
* `dsa_model=0`: Efficiency model from [Kang et. al. 2007][http://arxiv.org/abs/0704.1521v1]
* `dsa_model=1`: Efficiency model from [Kang & Ryu 2013][https://arxiv.org/abs/1212.3246]
* `dsa_model=2`: Efficiency model from [Ryu et. al. (2019)][https://arxiv.org/abs/1905.04476] extrapolated from the range `2.25 < M <= 5`. 
* `dsa_model=3`: Efficiency model from [Caprioli & Pitkovsky (2014)][https://iopscience.iop.org/article/10.1088/0004-637X/783/2/91] used as in [Vazza et. al. (2012)][https://arxiv.org/abs/1201.3362] as `0.5 * Kang&Ryu 2013`
* `dsa_model=4`: Constant efficiency from [Pfrommer et. al. (2016)][doi:10.1088/0004-637X/764/1/95]

### No pre-existing CRs in the shocktube 

A minimal working version for the solution of the CR shock discussed in Pfrommer+16 would be:

```julia
par = RiemannParameters(Pl=63.499, Pr=0.1, t=1.5, dsa_model=4)
```

This also returns a parameter object: [`SodCRParameters_noCRs`](@ref) which can be found in `cr_sod_shock_noprepopulation.jl`.

### Pre-existing CRs in the shocktube

!!! Under Construction !!!


## Solving the shock

To solve the shock with the given initial condition you just need to call

```julia
sol = solve(x, par)
```

with par being either of the above mentioned parameter objects, multiple dispatch will take care of the rest.

`x` has to be an array with either sample positions along the tube, or your actual particle positions, to make calculating errors easier. You can also just pass it an array with a single position, if you're only interested in that specific part of the shock ( e.g. `x = [86.0]` for the center of the postshock region.)

This will return a solution object depending on which shock you're solving.

For the pure hydro case this is:

```julia
mutable struct SodHydroSolution
    x::Array{Float64,1}         # array of given positions
    rho::Array{Float64,1}       # array of densities along the tube
    rho4::Float64               # density in postshock region
    rho3::Float64               # density between contact disc. and rarefaction wave
    P::Array{Float64,1}         # array of pressures along the tube
    P34::Float64                # pressure between shock and rarefaction wave
    U::Array{Float64,1}         # array of internal energies along the tube
    v::Array{Float64,1}         # array of velocities along the tube
    v34::Float64                # velocity between shock and rarefaction wave
    vt::Float64                 # velocity of rarefaction wave
    vs::Float64                 # shock velocity
    Mach::Float64               # Mach number
end
```

## Utility


A common issue is running into the error `DomainError` when solving a CR Sod shock.
This is due to the definition of the incomplete beta function. You can avoid this by supplying a value for `xs_first_guess`, which is a first guess for the value of the shock compression ratio.
In case you don't know the target `xs` (which is the usual case) and are tired of trying different values there's a helper function for that:

```julia
function find_xs_first_guess(Ul::Float64, Mach::Float64;
                             xs_start::Float64=3.8, delta_xs::Float64=1.e-4,
                             dsa_model::Int64=2, thetaB::Float64=0.0)

    [...]
end
```

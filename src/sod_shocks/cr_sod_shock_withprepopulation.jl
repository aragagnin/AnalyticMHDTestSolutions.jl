#using Roots
#using Distributions
using NLsolve
using QuadGK
#using ForwardDiff
using SpecialFunctions
#include("cr_sod_shock_main.jl")
#include("beta_inc_functions.jl")

mutable struct SodCRParameters_withCRs

    rhol::Float64
    rhor::Float64
    Pl::Float64
    Pr::Float64
    Ul::Float64
    Ur::Float64
    P_cr_l::Float64
    P_cr_r::Float64
    E_cr_l::Float64
    E_cr_r::Float64
    cl::Float64
    cr::Float64
    M::Float64
    t::Float64
    x_contact::Float64
    Pe_ratio::Float64
    γ_th::Float64
    γ_cr::Float64
    Δγ::Float64
    γ_exp::Float64
    α::Float64
    β::Float64
    η2::Float64
    ζ::Float64 
    ξ::Float64

    function SodCRParameters_withCRs(;rhol::Float64=1.0, rhor::Float64=0.125,
                                Pl::Float64=0.0,         Pr::Float64=0.0,
                                Ul::Float64=0.0,         Ur::Float64=0.0,
                                P_cr_l::Float64=0.0,     P_cr_r::Float64=0.0,
                                E_cr_l::Float64=0.0,     E_cr_r::Float64=0.0,
                                Mach::Float64=0.0,       t::Float64,
                                x_contact::Float64=70.0,
                                Pe_ratio::Float64=0.01,
                                thetaB::Float64=0.0,
                                theta_crit::Float64=0.25π,
                                dsa_model::Int64=-1,
                                γ_th::Float64=5.0/3.0,
                                γ_cr::Float64=4.0/3.0 )

        γ_exp    = ( γ_th - 1.0 )/( 2.0 * γ_th )
        η2       = ( γ_th - 1.0 )/( γ_th + 1.0 )
        Δγ       =   γ_th - γ_cr

        # calculate Ul and Pl depending on input
        if (Pl == 0.0) & (Ul != 0.0)
            Pl = ( γ_th - 1.0 ) * rhol * Ul
        elseif (Ul == 0.0) & (Pl != 0.0)
            Ul = Pl / ( (γ_th - 1.0) * rhol )
        else
            error("Both Ul and Pl are zero!")
        end

        # selection of DSA model structs
        if dsa_model == 0
            acc_model = KR07()
        elseif dsa_model == 1
            acc_model = KR13()
        elseif dsa_model == 2
            acc_model = Ryu19()
        elseif dsa_model == 3
            acc_model = CS14()
        elseif dsa_model == 4
            acc_model = P16()
        elseif dsa_model == 5
            acc_model = NullAcc()
        else
            error("Invalid DSA model selection!\n
                   Pick one of the available models, or solve a pure Hydro shock with:\n
                   SodParameters")
        end

        # calculate Ur and Pr depending on input
        if (Pr == 0.0) & (Ur != 0.0)
            Pr = ( γ_th - 1.0 ) * rhor * Ur
        elseif (Ur == 0.0) & (Pr != 0.0)
            Ur = Pr / ( (γ_th - 1.0) * rhor )
        elseif (Ur == 0.0) & (Pr == 0.0) & (Mach == 0.0)
            println("Error! Ur, Pr and Mach are zero! Can't find solution!")
        else
            println("Both Ur and Pr are zero! Will calculate them depending on Machnumber.")
            Pr = solvePrfromMach(rhol, rhor, Pl, Mach, γ_th, γ_cr, acc_model)
            Ur = Pr / ( (γ_th - 1.0) * rhor )
        end
        #
        # if Mach == 0.0
        #     Mach = solveMach(Pl, Pr, rhol, rhor, γ_th)
        # end


        # cl = √(γ_eff(Pl, P_cr_l, γ_th, γ_cr) * (Pl + P_cr_l) / rhol)
        # cr = √(γ_eff(Pr, P_cr_r, γ_th, γ_cr) * (Pr + P_cr_r) / rhor)

        cl = √( ( γ_th *  Pl + γ_cr * P_cr_l ) / rhol )
        cr = √( ( γ_th *  Pr + γ_cr * P_cr_r ) / rhor )

        α = (γ_cr - 1.0)/(2.0 * Δγ)
        β = (1.0 - γ_th)/(2.0 * Δγ)

        E_cr_r = P_cr_r/( ( γ_cr - 1.0 )*rhor )
        E_cr_l = P_cr_l/( ( γ_cr - 1.0 )*rhol )

        # calculate B angle dependent efficiency following Pais+ 2018, MNRAS, 478, 5278
        delta_theta = π/18.0
        thetaB *= (π/180.0)
        etaB = 0.5*( tanh( (theta_crit - thetaB)/delta_theta ) + 1.0 )

        Xcr = P_cr_r/Pr
        ζ = etaB*calc_η_Ms(acc_model, Mach, Xcr)

        ξ = ζ/(1.0 - ζ)

        new(rhol, rhor,
            Pl, Pr,
            Ul, Ur,
            P_cr_l, P_cr_r,
            E_cr_l, E_cr_r,
            cl, cr,
            Mach, t,
            x_contact,
            Pe_ratio,
            γ_th, γ_cr,
            Δγ,
            γ_exp,
            α, β,
            η2,
            ζ, ξ)
    end
end



"""
        Integral function
"""

"""
    A(P::Float64, ρ::Float64, γ::Float64

Entropy of gas (CP+16)
"""
@inline A(P::T, ρ::T, γ::Float64) where T = P * ρ^(-γ)

"""
    Ã(P::Float64, ρ::Float64, γ::Float64)

Entropy times adiabatic index (CP+16, after Eq. (C3) ).
"""
@inline Ã(P::T, ρ::T, γ::Float64) where T = γ * A(P, ρ, γ)

"""
    a(P::Float64, ρ::Float64, γ::Float64) 

Entropy times gamma times rho^gamma.
"""
@inline a(P::T, ρ::T, γ::Float64) where T = Ã(P, ρ, γ) * ρ^γ

"""
    x_func(ρ::Float64, P_th::Float64, P_cr::Float64, par::SodCRParameters_withCRs)

Base for incomplete gamma function (CP+16, Eq. C4 ).
"""
@inline function x_func(ρ::T, P_th::T, P_cr::T, par::SodCRParameters_withCRs) where T
    return  a(P_th, ρ, par.γ_th)/(a(P_th, ρ, par.γ_th) + a(P_cr, ρ, par.γ_cr))
end



@inline function reduced_Beta(x::T, α::T, β::T) where T
    return 1.0/beta(α, β) * x^(α - 1.0) * ( 1.0 - x )^(β - 1.0)
end

@inline function γ_eff(P_th::T, P_cr::T, γ_th::T, γ_cr::T) where T
    return (γ_cr*P_cr + γ_th*P_th)/(P_th + P_cr)
end

@inline function incomplete_beta(a, b, t) 
    return t^(a-1.0) * (1.0-t)^(b-1.0)
end

"""
    function I(ρ::T, P_th::T, P_cr::T, par::SodCRParameters_withCRs ) where T

Simplified Integral, Eq. C4 in CP+16
"""
@inline function I(ρ::T, P_th::T, P_cr::T, par::SodCRParameters_withCRs ) where T

    x_   = x_func(ρ, P_th, P_cr, par)
    Ã_cr = Ã(P_cr, ρ, par.γ_cr)
    Ã_th = Ã(P_th, ρ, par.γ_th)
    B(x) = incomplete_beta(par.α, par.β, x)

    result_B, result_error = quadgk(B, 0, x_, rtol=1.e-4)

    return √(Ã_cr)/par.Δγ * (Ã_cr/Ã_th)^par.α * result_B
end

# """
#     function I(ρ::T, P_th::T, P_cr::T, par::SodCRParameters_withCRs ) where T

# Simplified Integral, Eq. C4 in CP+16
# """
# @inline function I(ρ::T, P_th::T, P_cr::T, par::SodCRParameters_withCRs ) where T
#     x_   = x_func(ρ, P_th, P_cr, par)
#     ℬ   = reduced_Beta(x_, par.α, par.β)
#     Ã_cr = Ã(P_cr, ρ, par.γ_cr)
#     Ã_th = Ã(P_th, ρ, par.γ_th)

#     # B(x) = incomplete_beta(par.α, par.β, x)

#     # result_B, result_error = quadgk(B, 0, x_, rtol=1.e-4)


#     return √(Ã_cr)/par.Δγ * (Ã_cr/Ã_th)^par.α * ℬ
# end

# """
#     function I(ρ::T, P_th::T, P_cr::T, par::SodCRParameters_withCRs ) where T

# Simplified Integral, Eq. C4 in CP+16
# """
# @inline function I_integrated(ρ::T, P_th::T, P_cr::T, par::SodCRParameters_withCRs ) where T
#     #x_   = x_func(ρ, P_th, P_cr, par)
#     Ã_cr = Ã(P_cr, ρ, par.γ_cr)
#     Ã_th = Ã(P_th, ρ, par.γ_th)

#     I_in(x) = √( Ã_cr * x^(par.γ_cr-3) + Ã_th * x^(par.γ_th-3) )

#     result_B, result_error = quadgk(I_in, 0, ρ, rtol=1.e-4)


#     return result_B
# end



"""
    Pressure
"""

function solveP4(par::SodCRParameters_withCRs, sol::SodCRSolution)
    1.0 / ( 1.0 + par.ζ ) * ( sol.P34_tot + par.ξ * ( par.γ_cr - 1.0 ) / ( par.γ_th - 1.0 ) * par.Pr * sol.xs^par.γ_th - par.P_cr_r * sol.xs^par.γ_th )
end

function solveP3(par::SodCRParameters_withCRs, sol::SodCRSolution)
    # Solves the pressure of the middle state in a shocktube

    Pcr3 = par.P_cr_l * ( sol.rho3/par.rhol )^par.γ_cr
    P3 = par.Pl * ( sol.rho3/par.rhol )^par.γ_th

    return P3, Pcr3
end

function solveP2(x::Float64; par::SodCRParameters_withCRs)
    # solves the pressure along the rarefaction wave
    return par.Pl * ( -par.η2 * x / ( par.cl * par.t ) + ( 1.0 - par.η2 ) )^( 1.0/par.γ_exp )
end

function solveP(x::Float64, par::SodCRParameters_withCRs, sol::SodCRSolution)
    # returns P values according to position in shocktube

    if x <= -par.cl * par.t
        return (par.Pl + par.P_cr_l),
                par.Pl, par.P_cr_l,
                (1.0 - par.Pe_ratio)*par.P_cr_l,
                par.Pe_ratio*par.P_cr_l

    elseif -par.cl * par.t < x <= -sol.vt * par.t
        return 0.0, 0.0, 0.0, 0.0, 0.0

    elseif -sol.vt * par.t < x <= sol.v34 * par.t
        return (sol.P3_th + sol.P3_cr),
                sol.P3_th, sol.P3_cr,
                (1.0 - par.Pe_ratio)*sol.P3_cr,
                par.Pe_ratio*sol.P3_cr

    elseif sol.v34 * par.t < x <= sol.vs * par.t
        return (sol.P4_th + sol.P4_cr),
                sol.P4_th, sol.P4_cr,
                (1.0 - par.Pe_ratio)*sol.P4_cr,
                par.Pe_ratio*sol.P4_cr

    elseif sol.vs * par.t < x
        return (par.Pr + par.P_cr_r),
                par.Pr, par.P_cr_r,
                (1.0 - par.Pe_ratio)*par.P_cr_r,
                par.Pe_ratio*par.P_cr_r
    else
        error("Error! x = $x  t = $(par.t)")
        return 0.0, 0.0, 0.0, 0.0, 0.0
    end
end


"""
    Density
"""
function solveRho2(x::Float64, par::SodCRParameters_withCRs)
    # solves density along the rarefaction wave

    #f(rho) = I(rho) - I(par.rhol) + x/par.t + sqrt.( A())
    return 0.0

end

"""
    xs(ρ₄::T, ρᵣ::T) where T

Compression ratio over shock.
``x_s = \frac{ρ_4}{ρ_r}``
"""
xs(ρ₄::T, ρᵣ::T) where T = ρ₄/ρᵣ

"""
    xr(ρ₃::T, ρₗ::T) where T

Compression ratio over rarefaction wave.
``x_r = \\frac{ρ_3}{ρ_l}``
"""
xr(ρ₃::T, ρₗ::T) where T = ρ₃/ρₗ


function rho34_solver!(F, x, par::SodCRParameters_withCRs)

    

    xₛ  = xs(x[1], par.rhor)
    xᵣ  = xr(x[2], par.rhol)

    γ_inj = 4.0/3.0 # fix!

    # CP+16, Eq. C6
    P₃₄ = par.P_cr_l*xᵣ^par.γ_cr + par.Pl*xᵣ^par.γ_th

    # CP+16, Eq. C7
    ϵ_th_ad = par.Ur * xₛ^par.γ_th

    ϵ_cr4 = par.E_cr_r * xₛ^par.γ_cr

    # CP+16, Eq. C8
    ε₄ = ( (1.0 - par.ζ) * ( par.γ_th - 1.0 ) / ( γ_inj - 1.0 ) + par.ζ)^(-1)  * 
         ( P₃₄ / ( γ_inj - 1.0 ) + par.ξ * ϵ_th_ad - (par.γ_cr - 1.0) / ( γ_inj - 1.0 ) * ϵ_cr4 ) + 
         ϵ_cr4 - par.ξ * ϵ_th_ad

    #println("ε₄ = $ε₄")

    # total energy density right 
    ϵᵣ = par.Ur + par.E_cr_r

    P_tot_r = (par.Pr + par.P_cr_r)

    # CP+16, Eq. C5
    F[1] = ( P₃₄ - P_tot_r ) * ( xₛ - 1.0 ) - 
            par.rhor * xₛ * ( I(par.rhol, par.Pl, par.P_cr_l, par ) - 
                I(xᵣ*par.rhol, par.Pl, par.P_cr_l, par ) )^2

    F[2] = ( P₃₄ + P_tot_r ) * ( xₛ - 1.0 ) + 2*( xₛ*ϵᵣ - ε₄ )

    F
end


function solveRho34(par::SodCRParameters_withCRs)

    rho34_helper!(F, x) = rho34_solver!(F, x, par)

    #x_guess = 0.5*(par.rhol + par.rhor)
    initial_x = [0.5par.rhol, 0.5*par.rhol]
    nlsolve(rho34_helper!, initial_x, method = :newton)#, autodiff = :forward)

        #   rho3           rho4
    return initial_x[1], initial_x[2]
end


function solveRho(x::Float64, par::SodCRParameters_withCRs, sol::SodCRSolution)
    # returns P values according to position in shocktube

    if x <= -par.cl * par.t
        return par.rhol
    elseif -par.cl * par.t < x <= -sol.vt * par.t
        return solveRho2(x, par)
    elseif -sol.vt * par.t < x <= sol.v34 * par.t
        return sol.rho3
    elseif sol.v34 * par.t < x <= sol.vs * par.t
        return sol.rho4
    elseif sol.vs * par.t < x
        return par.rhor
    else
        error("Solution outside simulation domain!")
        return 0.0
    end
end


"""
    Velocity
"""
function solveV34(par::SodCRParameters_withCRs, sol::SodCRSolution)
    # solves velocity along the isopressure region 2-3
    return √( (sol.P34_tot - (par.Pr + par.P_cr_r ) ) * (sol.rho4 - par.rhor)/sol.rho4*par.rhor )
end

function solveV2(x::Float64, par::SodCRParameters_withCRs)
    # solves the velocity along the rarefaction wave
    return 0.0
end

function solveVs(par::SodCRParameters_withCRs, sol::SodCRSolution)
    # solves the shock velocity
    return sol.v34 * sol.rho4 / ( sol.rho4 - par.rhor )
end

function solveVt(par::SodCRParameters_withCRs, sol::SodCRSolution)
    # solves the velocity of the tail of the rarefaction wave
    return I(sol.rho3, sol.P3_th, sol.P3_cr, par) - I(par.rhol, par.Pl, par.P_cr_l, par) +
            sqrt(
                A(sol.P3_cr, sol.rho3, par.γ_cr)*sol.rho3^(par.γ_cr - 1.0) +
                A(sol.P3_th, sol.rho3, par.γ_th)*sol.rho3^(par.γ_th - 1.0) )
end

function solveV(x::Float64, par::SodCRParameters_withCRs, sol::SodCRSolution)
    # solves the velocity along the shocktube
    if x <= -par.cl * par.t
        return 0.0
    elseif -par.cl * par.t < x <= -sol.vt * par.t
        return solveV2(x, par)
    elseif -sol.vt * par.t < x <= sol.vs * par.t
        return sol.v34
    elseif sol.vs * par.t < x
        return 0.0
    else
        println("Error!")
        return 0.0
    end
end


"""
    Main function for solver
"""


function solveSodShockCR_withPrepopulation(x::Vector{Float64}; par::SodCRParameters_withCRs)

    # set up datatype to store riemann solution
    sol = SodCRSolution(x)

    println("par.P_cr_l = $(par.P_cr_l)")

    # transform into rest-frame of contact discontiuity
    x_in = sol.x .- par.x_contact

    # solve Pressure
    sol.rho3, sol.rho4 = solveRho34(par)

    sol.xs  = sol.rho4/par.rhor
    sol.xr  = sol.rho3/par.rhol 

    sol.P34_tot = par.P_cr_l*sol.xr^par.γ_cr + par.Pl*sol.xr^par.γ_th

    sol.P3_cr = par.P_cr_l*sol.xr^par.γ_cr
    sol.P3_th = sol.P34_tot - sol.P3_cr

    # P_inj     = P_inj_f(sol.xs, par.rhor, par.Pr, par.γ_th, par.γ_cr, par.ξ)
    P_inj = 0.0
    # sol.P4_cr = par.P_cr_r*sol.xs^par.γ_cr + P_inj
    
    sol.P4_th = solveP4(par, sol)
    sol.P4_cr = sol.P34_tot - sol.P4_th


    println("P_inj = $P_inj\tsol.P4_cr = $(sol.P4_cr)\tsol.P4_th = $(sol.P4_th)")

    # solve velocity
    sol.v34 = solveV34(par, sol)

    sol.vs = solveVs(par, sol)
    sol.vt = solveVt(par, sol)

    for i=1:length(sol.x)
        sol.v[i] = solveV(x_in[i], par, sol)
    end

    for i = 1:length(sol.x)
        sol.P_tot[i], sol.P_th[i], sol.P_cr[i], sol.P_cr_p[i], sol.P_cr_e[i] = solveP(x_in[i], par, sol)
    end

    for i = 1:length(sol.x)
        sol.rho[i] = solveRho(x_in[i], par, sol)
    end

    sol.Mach = sol.vs/par.cr

    sol.U_tot = sol.P_tot ./ ( (par.γ_th - 1.0) .* sol.rho )
    sol.U_th = sol.P_th ./ ( (par.γ_th - 1.0) .* sol.rho )
    sol.E_cr_p = sol.P_cr_p ./ ( (par.γ_cr - 1.0) .* sol.rho )
    sol.E_cr_e = sol.P_cr_e ./ ( (par.γ_cr - 1.0) .* sol.rho )

    println("par.P_cr_r = $(par.P_cr_r)")

    return sol
end

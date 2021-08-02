"""
    Types for acceleration efficiencies
"""

"""
    ShockAccelerationEfficiency

Abstract type for shock acceleration efficiencies
"""
abstract type ShockAccelerationEfficiency end

"""
    struct KR07{T} <: ShockAccelerationEfficiency
        X::T
        η_max::T
    end
"""
struct KR07{T} <: ShockAccelerationEfficiency
    X::T
    η_max::T

    """
        KR07(X::T=0.05, η_max::T=0.0348) where T

    Default values for Kang, Ryu, Cen, Ostriker 2007, http://arxiv.org/abs/0704.1521v1
    """
    KR07(X::T=0.3, η_max::T=0.57) where T = new{T}(X, η_max)
end




"""
    struct KR13{T} <: ShockAccelerationEfficiency
        X::T
        η_max::T
    end
"""
struct KR13{T} <: ShockAccelerationEfficiency
    X::T
    η_max::T

    """
        KR13(X::T=0.05, η_max::T=0.0348) where T

    Default values for Kang&Ryu 2013 efficiancy.
    """
    KR13(X::T=0.05, η_max::T=0.2055) where T = new{T}(X, η_max)
end




"""
    struct Ryu19{T} <: ShockAccelerationEfficiency
        X::T
        η_max::T
    end
"""
struct Ryu19{T} <: ShockAccelerationEfficiency
    X::T
    η_max::T

    """
        Ryu19(X::T=0.05, η_max::T=0.0348) where T 

    Default values for Ryu+19 efficiancy.
    """
    Ryu19(X::T=0.05, η_max::T=0.0348) where T = new{T}(X, η_max)
end




"""
    struct CS14{T} <: ShockAccelerationEfficiency
        X::T
        η_max::T
    end
"""
struct CS14{T} <: ShockAccelerationEfficiency
    X::T
    η_max::T

    """
        CS14(X::T=0.05, η_max::T=0.5*0.2055) where T

    Default values for Caprioli&Spitkovsky 2014 efficiancy.
    """
    CS14(X::T=0.05, η_max::T=0.5*0.2055) where T = new{T}(X, η_max)
end




"""
    struct P16{T} <: ShockAccelerationEfficiency
        X::T
        η_max::T
    end
"""
struct P16{T} <: ShockAccelerationEfficiency
    X::T
    η_max::T

    """
        P16(X::T=0.05, η_max::T=0.5) where T

    Default values for constant efficiency as in Pfrommer+ 2016, doi: 10.1093/mnras/stw2941 
    """
    P16(X::T=0.05, η_max::T=0.5) where T = new{T}(X, η_max)
end

"""
    struct P16{T} <: ShockAccelerationEfficiency
        X::T
        η_max::T
    end
"""
struct NullAcc{T} <: ShockAccelerationEfficiency
    X::T
    η_max::T

    """
        NullAcc(X::T=0.05, η_max::T=0.5) where T

    Fallback option with no injection.
    """
    NullAcc(X::T=0.0, η_max::T=0.0) where T = new{T}(X, η_max)
end


"""
    Functions
"""


"""
    Helper functions
"""
# power functions
@inline p2(x::T) where T = x*x
@inline p3(x::T) where T = x*x*x

"""
    kr_fitting_function(x::Real, 
                             a0::Real, a1::Real, a2::Real, a3::Real, a4::Real)

Helper function to use the fitting function from KR07.
"""
@inline function kr_fitting_function(x::Real, 
                                    a0::Real, a1::Real, a2::Real, a3::Real, a4::Real)
    mm = x - 1.0

	return ( a0 + a1*mm + a2*p2(mm) + a3*p3(mm) + a4*p2(p2(mm)) ) / p2(p2(x))
end



"""
    KR07_acc(M::T) where T

Efficiency model from Kang, Ryu, Cen, Ostriker 2007, http://arxiv.org/abs/0704.1521v1
"""
function η_Ms_acc(efficiancy::KR07, M::T) where T
    if M < 1.0
        return 0.0
    elseif M <= 2.0
        return 1.96e-3*(M^2 - 1.)             # eq. A3
    else
        return kr_fitting_function(M, 5.46, -9.78, 4.17, -0.334, 0.57)
    end
end

"""
    KR07_reacc(M::T) where T

Efficiency model from Kang, Ryu, Cen, Ostriker 2007, http://arxiv.org/abs/0704.1521v1
"""
function η_Ms_reacc(efficiancy::KR07, M::T) where T
    if M <= 1.5
        # adiabatic index of gas
        γ::T  = 5.0/3.0
        m2::T = M^2
        # analytic compression ratio
        xs::T = (γ + 1.0) / ( γ -1.0 + 0.5m2 )
        δ0::T = 2.0 * ( (2γ*m2 - γ + 1.0) / (γ + 1.0) - (xs^γ) ) / 
                ( γ * (γ - 1.0) * m2 * xs )
        return 1.025*δ0
    else
        return kr_fitting_function(M, 0.24, -1.56, 2.8, 0.512, 0.557)
    end
end



"""
    η_Ms_acc(efficiancy::KR13, M::T) where T

Efficiency model from Kang&Ryu 2013, doi:10.1088/0004-637X/764/1/95
"""
function η_Ms_acc(efficiancy::KR13, M::T) where T

    if M < 2.0
        return 0.0
    elseif 2.0 <= M <= 5.0
        return -0.0005950569221922047 + 1.880258286365841e-5 * M^5.334076006529829 
    elseif 5.0 < M <= 15.0
        return kr_fitting_function(M, -2.8696966498579606, 9.667563166507879,
                                  -8.877138312318019, 1.938386688261113, 0.1806112438315771)
    else
        return 0.21152
    end
end

"""
    η_Ms_reacc(efficiancy::KR13, M::T) where T

Efficiency model from Kang&Ryu 2013, doi:10.1088/0004-637X/764/1/95
"""
function η_Ms_reacc(efficiancy::KR13, M::T) where T

    if M < 2.0
        return 0.0
    elseif M <= 17.7
        return kr_fitting_function(M, -0.722, 2.7307, -3.2854, 1.3428, 0.1901)
    else
        return 0.2055
    end
end


"""
    η_Ms_acc(efficiancy::Ryu19, M::T) where T

Efficiency model from Ryu et al. 2019, https://arxiv.org/abs/1905.04476
Values for 2.25 < M <= 5.0 extrapolated to entire range
"""
function η_Ms_acc(efficiancy::Ryu19, M::T) where T

    if M < 2.25
        return 0.0
    elseif M <= 34.0
        return kr_fitting_function(M, -1.5255, 2.4026, -1.2534, 0.2215, 0.0336)
    else
        return 0.0348
    end
end


"""
    η_Ms_reacc(efficiancy::Ryu19, M::T) where T

Efficiency model from Ryu et al. 2019, https://arxiv.org/abs/1905.04476
Values for 2.25 < M <= 5.0 extrapolated to entire range
"""
function η_Ms_reacc(efficiancy::Ryu19, M::T) where T

    if M < 2.25
        return 0.0
    elseif M <= 34.0
        return kr_fitting_function(M, 0.3965, -0.21898, -0.2074, 0.1319, 0.0351)
    else
        return 0.0348
    end
end


"""
    η_Ms_acc(efficiancy::CS14, M::T) where T

Efficiency from Caprioli&Spitkovsky 2015, doi: 10.1088/0004-637x/783/2/91
Same simplified approach as Vazza+12 -> is roughly half the efficiency of Kang&Ryu 2013.
"""
function η_Ms_acc(efficiancy::CS14, M::T) where T
    vazza_factor = 0.5
    return vazza_factor * η_Ms_acc(KR13(), M)
end

"""
    η_Ms_reacc(efficiancy::CS14, M::T) where T

Efficiency from Caprioli&Spitkovsky 2015, doi: 10.1088/0004-637x/783/2/91
Same simplified approach as Vazza+12 -> is roughly half the efficiency of Kang&Ryu 2013.
"""
function η_Ms_reacc(efficiancy::CS14, M::T) where T
    vazza_factor = 0.5
    return vazza_factor * η_Ms_reacc(KR13(), M)
end

"""
    P16_acc(M::T) where T

Constant efficiency as in Pfrommer+ 2016, doi: 10.1093/mnras/stw2941 
"""
function P16_acc(M::T) where T
    return 0.5
end

"""
    η_Ms_acc(efficiancy::P16, M::T) where T

Constant efficiency as in Pfrommer+ 2016, doi: 10.1093/mnras/stw2941 
"""
function η_Ms_acc(efficiancy::P16, M::T) where T
    return 0.5
end

"""
    η_Ms_reacc(efficiancy::P16, M::T) where T

Constant efficiency as in Pfrommer+ 2016, doi: 10.1093/mnras/stw2941 
"""
function η_Ms_reacc(efficiancy::P16, M::T) where T
    return 0.5
end

"""
    η_Ms_acc(efficiancy::NullAcc, M::T) where T

Fallback with no injection.
"""
function η_Ms_acc(efficiancy::NullAcc, M::T) where T
    return 0.0
end

"""
    η_Ms_reacc(efficiancy::NullAcc, M::T) where T

Fallback with no injection.
"""
function η_Ms_reacc(efficiancy::NullAcc, M::T) where T
    return 0.0
end


"""
    calc_η_Ms(η_model::ShockAccelerationEfficiency, M::T, X::T) where T

Calculates the efficiancy as a linear interpolation between acceleration and reacceleration.
"""
function calc_η_Ms(η_model::ShockAccelerationEfficiency, M::T, X::T) where T

    # check if pressure ratio is greater than the target ratio
    if X > η_model.X 
        X = η_model.X
    end

    # initial acceleration
    η1 = η_Ms_acc(η_model, M)
    X1::T = 0.0

    # reacceleration
    η2 = η_Ms_reacc(η_model, M)
    X2 = η_model.X

    # linear interpolation between the two cases:
    η_tot::T = η1 + (X - X1) * (η2 - η1) / (X2 - X1)

    # check if η_tot is smaller than the maximum efficiancy
    if η_tot > η_model.η_max
        return η_model.η_max
    else
        return η_tot
    end

end


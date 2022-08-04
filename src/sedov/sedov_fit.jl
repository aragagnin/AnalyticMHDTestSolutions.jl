using SpecialFunctions: gamma
using Interpolations: linear_interpolation

struct w1Data
    ndim::Int64
    γ::Float64

    ω::Vector{Float64}
    β::Vector{Float64}
    C::Vector{Float64}
end

struct w2Data
    ndim::Int64
    γ::Float64
    
    ω::Vector{Float64}
    β::Vector{Float64}
    C::Vector{Float64}
end

struct w3Data
    ndim::Int64
    γ::Float64
    
    ω::Vector{Float64}
    β::Vector{Float64}
    C::Vector{Float64}
end

struct elseData
    ndim::Int64
    γ::Float64
    
    ω::Vector{Float64}
    β::Vector{Float64}
    C::Vector{Float64}
end

"""
        w1 functions
"""
function xifunc(F::Float64, c::w1Data)
    return F^(-c.β[7]) * (c.C[2] * (c.F - c.C[3]))^c.β[3] * (c.C[4] * (c.C[5] - F))^(-c.β[2])
end
function Dfunc(F::Float64, c::w1Data)
    return xifunc(F, c)^(c.ndim-2)
end
function Vfunc(F::Float64, c::w1Data)
    return xifunc(F, c)
end
function Pfunc(F::Float64, c::w1Data)
    return xifunc(F, c)^c.ndim
end

"""
        w2 functions
"""
function xifunc(F::Float64, c::w2Data)
    return  F^(-c.β[7]) * (c.C[2] * (F - c.C[3]))^((c.γ - 1.0 )*c.β[1]) *
            exp(((c.γ + 1.0)*c.β[1]*(1.0 - F)) / (F - c.C[3]))
end
function Dfunc(F::Float64, c::w2Data)
    return  F^c.β[8] * (c.C[2] * (F - c.C[3]))^((4.0 - c.ndim - 2.0*c.γ)*c.β[]) *
            (c.C[6] * (c.C[7] - F))^(-c.β[6]) *
            exp((-2.0*(c.γ + 1.0)*c.β[1]*(1.0 - F)) / (F - c.C[3]))
end
function Vfunc(F::Float64, c::w2Data)
    return xifunc(F, c) * F
end
function Pfunc(F::Float64, c::w2Data)
    return F^c.b8 * (c.C[2] * (F - c.C[3]))^(-c.ndim * c.γ * c.β[1]) *
           (c.C[6] * (c.C[7] - F))^(1.0 - c.β[6])
end

"""
        w3 functions
"""
function xifunc(F::Float64, c::w3Data)
    return F^(-c.β[7]) * (c.C[2] * (F - c.C[3]))^c.β[3] * (c.C[6] * (c.C[7] - F))^(-c.β[2])
end
function Dfunc(F::Float64, c::w3Data)
    return F^c.β[8] * (c.C[2] * (F - c.C[3]))^(c.β[4] - c.ω[1] * c.β[3]) *
           (c.C[6] * (c.C[7] - F))^(1.0 - 4.0*c.β[1]) *
           exp((-c.ndim * c.γ * (c.γ + 1.0) * c.β[1] * (1.0 - F)) /
               (c.C[7] - F))
end
function Vfunc(F::Float64, c::w3Data)
    return xifunc(F, c) * F
end
function Pfunc(F::Float64, c::w3Data)
    return F^c.b8 * (c.C[6] * (c.C[7] - F))^(2.0*(c.ndim * c.γ - c.ndim - c.γ )*c.β[1]) *
           exp((-c.ndim * c.γ * ( c.γ + 1.0) * c.β[1] * (1.0 - F)) / (c.C[7] - F))
end

"""
    else functions
"""
function xifunc(F::Float64, c::elseData)
    return F^(-c.β[7]) * (c.C[2] * (F - c.C[3]))^c.β[3] * (c.C[4] * (c.C[5] - F))^(-c.β[2])
end
function Dfunc(F::Float64, c::elseData)
    return F^c.β[8] * (c.C[2] * (F - c.C[3]))^(c.β[4] - c.ω[1] * c.β[3]) *
           (c.C[4] * (c.C[5] - F))^(c.β[5] + c.ω[1] * c.β[2]) *
           (c.C[6] * (c.C[7] - F))^(-c.β[6])
end
function Vfunc(F::Float64, c::elseData)
    return xifunc(F, c) * F
end
function Pfunc(F::Float64, c::elseData)
    return F^(c.β[9]) * (c.C[4] * (c.C[5] - F))^(c.β[5] + (c.ω[1] - 2.0)*c.β[2]) *
           (c.C[6] * (c.C[7] - F))^(1.0 - c.β[6])
end



function get_ideal_sedov_fit(ndim::Integer=3; γ::Real=5.0/3.0)

    ω = zeros(4)
    
    # ω[1] parametrizes the background density as rho = A * r^(-ω[1])
    ω[1] = 0.0 # constant background density

    # The following is copied from Book (1991)
    # The omega parameters
    ω[2] = (3.0 * ndim - 2.0 + γ*( 2.0 - ndim )) / (γ + 1.0)
    ω[3] = (2.0 * ( γ - 1.0 ) + ndim)  / γ
    ω[4] = ndim * (2.0 - γ)

    # The beta_n parameters:
    β = zeros(9)
    β[1] = 1.0 / (ndim * γ - ndim + 2)
    β[3] = (γ - 1.0) / ( γ * ( ω[3] - ω[1] ))
    β[4] = (ndim - ω[1]) / (γ*(ω[3]-ω[1]))
    β[6] = ( 2.0*ndim - ω[1] * ( γ + 1 )) / ( ω[4] - ω[1] )
    β[7] = 2.0 / (ndim + 2.0 - ω[1])
    β[2] = β[3] + ( γ + 1 ) * β[1] - β[7]
    β[5] = β[2] * ((ndim - ω[1]) * (ndim + 2 - ω[1])) / (ω[4] - ω[1])
    β[8] = ω[1] * β[7]
    β[9] = ndim * β[7]
    
    # The C parameters
    C = zeros(7)
    C[1] = 2^ndim * π^((ndim - 1.0) / 2.0 ) * gamma((ndim+1)/2) / gamma(ndim)
    C[6] = 2.0 / ( γ - 1.0 )
    C[7] = ( γ + 1.0 ) * 0.5
    C[2] = γ * C[6]
    C[3] = C[7] / γ
    C[4] = ( ndim * γ - ndim + 2.0 ) / ( (ω[2] - ω[1])*C[7])
    C[5] = ( ndim + 2.0 - ω[1]) * β[1] * C[7]
    # Eq. (8) to (11) in Book (1991) or their later alternatives

    if ω[1] == ω[2]
        c = w1Data(ndim, γ,
                   ω, β, C)
    elseif ω[1] == ω[3]
        c = w2Data(ndim, γ,
                   ω, β, C)
    elseif ω[1] == ω[4]
        c = w3Data(ndim, γ,
                   ω, β, C)
    else
        c = elseData(ndim, γ,
                     ω, β, C)
    end


    if ω[2] > ω[1]
        Fmin = C[3]
    else
        Fmin = C[7]
    end

    # We use F as a proxy to do the integration over xi
    F = range(Fmin, stop=1.0, length=100_000)
    xi = zeros(length(F))
    for i = 1:length(F)
        xi[i] = xifunc(F[i], c)
    end

    k = sortperm(F)
    F = F[k] # sort F according to the xi-values it produces
    xi = xi[k]

    D = zeros(length(F))
    V = zeros(length(F))
    P = zeros(length(F))
    for i = 1:length(F)
        D[i] = Dfunc(F[i], c)
        V[i] = Vfunc(F[i], c)
        P[i] = Pfunc(F[i], c)
    end

    D_interpolated = linear_interpolation(xi, D)
    P_interpolated = linear_interpolation(xi, P)
    V_interpolated = linear_interpolation(xi, V)
    I = @. xi^(ndim-1) * (D*V^2 + P) # Integrand
    dxi = zeros(length(I))
    # Use middle points
    for i = 1:(length(I)-1)
        I[i] = 0.5 * (I[i] + I[i+1])
        dxi[i] = xi[i+1] - xi[i]
    end
    integral = sum( @. I * dxi )
    # Finally, get alpha
    alpha = ((8.0 * C[1]) / ((γ^2 - 1.0)*(ndim + 2.0 - ω[1])^2)) * integral

    return SedovFit(alpha, ω[1], D_interpolated, P_interpolated, V_interpolated)
end
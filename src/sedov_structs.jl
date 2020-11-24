struct SedovData

     t::Float64
     m::Float64
     r::Vector{Float32}
     vr::Vector{Float32}
     rho::Vector{Float32}
     U::Vector{Float32}
     P::Vector{Float32}
     Pth::Vector{Float32}
     Pcr::Vector{Float32}
     E::Float64
     hsml::Vector{Float32}
     mach::Vector{Float32}
     r_shock_rho::Float64
     r_shock_P::Float64
     rho_out::Float64
     rho_s::Float64
     cs_out::Float64
     γ::Float64
     ndim::Int64

     function SedovData(t::Float64,             m::Float64,
                        r::Vector{Float32},     vr::Vector{Float32},
                        rho::Vector{Float32},   U::Vector{Float32},
                        Pcr::Vector{Float32},   E::Float64,
                        hsml::Vector{Float32},  mach::Vector{Float32};
                        γ::Float64=5.0/3.0,     ndim::Int64=3)

        P = @. (γ - 1.0 ) * rho * U
        Pth = P

        if Pcr != zeros(Float32, length(P))
            P += Pcr
        end


        firstbin = Int(floor(length(rho) * 0.99)-1)

        rho_out = mean( rho[ firstbin:end ])
        P_out   = mean( P[ firstbin:end ] )

        cs_out  = sqrt( 5.0/3.0 * P_out/rho_out)

        rho_s = rho_out * ( γ + 1.0 )/( γ - 1.0 )

        k = findmax(rho)[2]
        r_shock_rho = r[k]

        k = findmax(P)[2]
        r_shock_P   = r[k]

        new(t, m, r, vr, rho, U, P, Pth, Pcr, E, hsml, mach, r_shock_rho, r_shock_P,
            rho_out, rho_s, cs_out, γ, ndim)
    end
end

struct SedovParameters
    ndim::Int64
    t::Float64
    E::Float64
    rho_out::Float64
    rho_s::Float64
    cs_out::Float64
    γ::Float64

    function SedovParameters(time::Real, Etot::Real=1.0, 
                             U_background::Real=0.01, rho_background::Real=1.0;
                             Ndim::Integer=3, γ::Real=5.0/3.0)

        rho_s  = rho_background * ( γ + 1.0 )/( γ - 1.0 )
        P_out  = ( γ - 1.0 ) * rho_background * U_background
        cs_out = √( γ * P_out / rho_background )

        new(Int64(Ndim), Float64(time), Float64(Etot), Float64(rho_background),
            Float64(rho_s), Float64(cs_out), Float64(γ))
    end
end

struct SedovFit

    alpha::Float64
    w::Float64
    D
    P
    V

end

struct SedovSolution

    r::Vector{Float32}
    vr::Vector{Float64}
    rho::Vector{Float64}
    P::Vector{Float64}
    U::Vector{Float64}
    r_shock::Float64
    v_shock::Float64
    xs::Float64
    mach::Float64
    alpha::Float64

end
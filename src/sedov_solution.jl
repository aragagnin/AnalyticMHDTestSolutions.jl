using GadgetIO
using Statistics
using ProgressMeter


function R_s_analytic(data::SedovParameters, fit::SedovFit)
    return (data.E * data.t^2 / ( fit.alpha * data.rho_out ))^(1.0/ ( 2.0 + data.ndim ))
end

function R_s_dot_analytic(data::SedovParameters, fit::SedovFit)
    t = maximum([data.t, 1.e-30])
    return 2.0 * R_s_analytic(data, fit) / ( (data.ndim + 2.0 - fit.w ) * t)
end

function v_s_analytic(data::SedovParameters, fit::SedovFit)
    return 2.0 * R_s_dot_analytic(data, fit) /
            ( data.γ + 1.0 )
end

function P_s_analytic( data::SedovParameters, fit::SedovFit )
    return 2.0 * data.rho_out * R_s_dot_analytic( data, fit)^2 / ( data.γ + 1.0 )
end

function rho_analytic(xi::Real, data::SedovParameters, fit::SedovFit)
    if xi < 1.0
        return data.rho_s * fit.D(xi)
    else
        return data.rho_out
    end
end

function P_analytic(xi::Real, Pₛ::Real, fit::SedovFit)

    if xi < 1.0
        return Pₛ * fit.P(xi)
    else
        return 0.0
    end
end

function v_analytic(xi::Float64, vₛ::Real, fit::SedovFit)
    if xi < 1.0
        return vₛ * fit.V(xi)
    else
        return 0.0
    end
end




function get_sedov_solution(r::Vector{<:Real}, par::SedovParameters, fit::SedovFit)

    N   = length(r)
    rho = zeros(N)
    P   = zeros(N)
    v   = zeros(N)

    rₛ = R_s_analytic(par, fit)
    vₛ = v_s_analytic(par, fit)
    Pₛ = P_s_analytic(par, fit)

    for i = 1:N
        xi = r[i]/rₛ
        rho[i] = rho_analytic(xi, par, fit)
        P[i]   = P_analytic(xi, Pₛ, fit)
        v[i]   = v_analytic(xi, vₛ, fit)
    end

    U = @. P / ( (par.γ - 1.0 ) * rho )

    xs = par.rho_s/par.rho_out

    mach = vₛ/par.cs_out

    return SedovSolution(r, v, rho, P, U, rₛ, vₛ, xs, mach, fit.alpha)
end

function solve(r::Vector{<:Real}, sedov_par::SedovParameters)


    sedov_fit  = get_ideal_sedov_fit( sedov_par.ndim, γ=sedov_par.γ )

    sedov_ideal = get_sedov_solution(r, sedov_par, sedov_fit)

    return sedov_ideal

end

function get_sedov_data_from_gadget(fi::String, blast_center::Vector{Float64}=[3.0, 3.0, 3.0];
                                      CRs=false, Nbins::Int64=0)

    h = head_to_obj(fi)

    t = h.time
    m = h.massarr[1]

    x = read_block(fi, "POS", parttype=0)

    r = @. sqrt( (x[:,1] - blast_center[1])^2 +
                 (x[:,2] - blast_center[2])^2 +
                 (x[:,3] - blast_center[3])^2 )

    k = sortperm(r)

    r = Float32.(r[k])

    v = read_block(fi, "VEL", parttype=0)[k,:]

    vr = @. sqrt( v[:,1]^2 + v[:,2]^2 + v[:,3]^2 )

    rho = read_block(fi, "RHO", parttype=0)[k, 1]

    U = read_block(fi, "U", parttype=0)[k, 1]

    hsml = read_block(fi, "HSML", parttype=0)[k, 1]

    try
        mach = read_block(fi, "MACH", parttype=0)[k, 1]
    catch
        mach = zeros(Float32, length(U))
    end

    if CRs
        CRpP = read_block(fi, "CRpP", parttype=0)[k, 1]
        Ecr = @. CRpP/(1.0/3.0 * rho)
        γ = 7.0/5.0
    else
        CRpP = zeros(Float32, length(U))
        Ecr = zeros(Float32, length(U))
        γ = 5.0/3.0
    end

    E = sum( @. (Ecr * m + U * m + 0.5 * m * vr^2) )

    if Nbins != 0
        N = length(r)

        rbinwidth = maximum(r)/Nbins
        r_raw    = r
        v_raw    = vr
        rho_raw  = rho
        U_raw    = U
        CRpP_raw = CRpP
        hsml_raw  = hsml
        mach_raw  = mach

        r     = zeros(Float32, Nbins+1)
        vr    = zeros(Float32, Nbins+1)
        rho   = zeros(Float32, Nbins+1)
        U     = zeros(Float32, Nbins+1)
        CRpP  = zeros(Float32, Nbins+1)
        hsml  = zeros(Float32, Nbins+1)
        mach  = zeros(Float32, Nbins+1)
        Npart = zeros(Float32, Nbins+1)

        @showprogress "Binning..." for i = 1:N
            bin = floor(Int64, r_raw[i]/rbinwidth) + 1

            vr[bin]     += v_raw[i]
            rho[bin]    += rho_raw[i]
            U[bin]      += U_raw[i]
            CRpP[bin]   += CRpP_raw[i]
            hsml[bin]  += hsml_raw[i]
            mach[bin]  += mach_raw[i]
            Npart[bin]  += 1
        end

        for i = 1:(Nbins+1)
            r[i] = (i-1) * rbinwidth
            vr[i]   /= Npart[i]
            rho[i]  /= Npart[i]
            U[i]    /= Npart[i]
            CRpP[i] /= Npart[i]
            mach[i] /= Npart[i]
            hsml[i] /= Npart[i]
        end

    end

    k = findall(isnan.(rho))

    deleteat!(r, k)
    deleteat!(vr, k)
    deleteat!(rho, k)
    deleteat!(U, k)
    deleteat!(CRpP, k)
    deleteat!(hsml, k)
    deleteat!(mach, k)


    return SedovData(t, m, r, vr, rho, U, CRpP, E, hsml, mach, γ=γ)

end

function get_sedov_data_from_arepo(fi::String, blast_center::Vector{Float64}=[3.0, 3.0, 3.0];
                                      CRs=false, Nbins::Int64=0)

    h = head_to_obj(fi)

    t = h.time
    m = h.massarr[1]

    x = read_block_by_name(fi, "POS", parttype=0)

    r = @. sqrt( (x[:,1] - blast_center[1])^2 +
                 (x[:,2] - blast_center[2])^2 +
                 (x[:,3] - blast_center[3])^2 )

    k = sortperm(r)

    r = Float32.(r[k])

    v = read_block_by_name(fi, "VEL", parttype=0)[k,:]

    vr = @. sqrt( v[:,1]^2 + v[:,2]^2 + v[:,3]^2 )

    rho = read_block_by_name(fi, "RHO", parttype=0)[k, 1]

    U = read_block_by_name(fi, "U", parttype=0)[k, 1]

    hsml = read_block_by_name(fi, "HSML", parttype=0)[k, 1]

    try
        mach = read_block_by_name(fi, "MACH", parttype=0)[k, 1]
    catch
        mach = zeros(Float32, length(U))
    end

    if CRs
        CRpP = read_block_by_name(fi, "CRpP", parttype=0)[k, 1]
        Ecr = @. CRpP/(1.0/3.0 * rho)
        γ = 7.0/5.0
    else
        CRpP = zeros(Float32, length(U))
        Ecr = zeros(Float32, length(U))
        γ = 5.0/3.0
    end

    E = sum( @. (Ecr * m + U * m + 0.5 * m * vr^2) )

    if Nbins != 0
        N = length(r)

        rbinwidth = maximum(r)/Nbins
        r_raw    = r
        v_raw    = vr
        rho_raw  = rho
        U_raw    = U
        CRpP_raw = CRpP
        hsml_raw  = hsml
        mach_raw  = mach

        r     = zeros(Float32, Nbins+1)
        vr    = zeros(Float32, Nbins+1)
        rho   = zeros(Float32, Nbins+1)
        U     = zeros(Float32, Nbins+1)
        CRpP  = zeros(Float32, Nbins+1)
        hsml  = zeros(Float32, Nbins+1)
        mach  = zeros(Float32, Nbins+1)
        Npart = zeros(Float32, Nbins+1)

        @showprogress "Binning..." for i = 1:N
            bin = floor(Int64, r_raw[i]/rbinwidth) + 1

            vr[bin]     += v_raw[i]
            rho[bin]    += rho_raw[i]
            U[bin]      += U_raw[i]
            CRpP[bin]   += CRpP_raw[i]
            hsml[bin]  += hsml_raw[i]
            mach[bin]  += mach_raw[i]
            Npart[bin]  += 1
        end

        for i = 1:(Nbins+1)
            r[i] = (i-1) * rbinwidth
            vr[i]   /= Npart[i]
            rho[i]  /= Npart[i]
            U[i]    /= Npart[i]
            CRpP[i] /= Npart[i]
            mach[i] /= Npart[i]
            hsml[i] /= Npart[i]
        end

    end

    k = findall(isnan.(rho))

    deleteat!(r, k)
    deleteat!(vr, k)
    deleteat!(rho, k)
    deleteat!(U, k)
    deleteat!(CRpP, k)
    deleteat!(hsml, k)
    deleteat!(mach, k)


    return SedovData(t, m, r, vr, rho, U, CRpP, E, hsml, mach, γ=γ)

end

function get_sedov_solution_from_gadget(filename::String, blast_center::Vector{Float64}=[3.0, 3.0, 3.0];
                            CRs::Bool=false, Nbins::Int64=0, Ndim::Int64=3)

    sedov_data = get_sedov_data_from_gadget(filename, blast_center, CRs=CRs, Nbins=Nbins)

    if CRs
        γ=7.0/5.0
    else
        γ=5.0/3.0
    end

    firstbin = Int(floor(length(sedov_data.rho) * 0.99)-1)

    U_out   = mean( sedov_data.U[ firstbin:end ] )

    sedov_par = SedovParameters(sedov_data.t, sedov_data.E, 
                             U_out, sedov_data.rho_out,
                             Ndim=Ndim, γ=γ)

    sedov_ideal = solve(sedov_data.r, sedov_par)

    return sedov_data, sedov_ideal

end


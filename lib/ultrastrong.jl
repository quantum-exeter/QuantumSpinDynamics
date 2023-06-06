########################
#### ultrastrong.jl ####
########################

function szTzero(prm::Lorentzian, ang::CouplingAngles, n::Levels)
    H = HTot(prm, ang, n)
    F = eigen(H)
    ϵ, P = F.values, F.vectors
    E_ground = ϵ[1]
    sz = 0
    σz′ = kronecker(σz, 𝕀(n.n1), 𝕀(n.n2))
    for i in 1:length(ϵ)
        if (ϵ[i] - E_ground) < 1e-16
            sz += tr(σz′*P[:,i]*adjoint(P[:,i]))
        end
    end
    return sz
end

function szTzero′(prm::Lorentzian, ang::CouplingAngles, n::Levels)
    H = HTot′(prm, ang, n)
    F = eigen(H)
    ϵ, P = F.values, F.vectors
    E_ground = ϵ[1]
    sz = 0
    σz′ = kronecker(σz, 𝕀(n.n1), 𝕀(n.n2))
    for i in 1:length(ϵ)
        if (ϵ[i] - E_ground) < 1e-16
            sz += tr(σz′*P[:,i]*adjoint(P[:,i]))
        end
    end
    return sz
end

function nTzero′(prm::Lorentzian, ang::CouplingAngles, n::Levels)
    H = HTot′(prm, ang, n)
    F = eigen(H)
    ϵ, P = F.values, F.vectors
    E_ground = ϵ[1]
    n1 = 0
    n2 = 0
    N1′ = kronecker(𝕀(2), N(n.n1), 𝕀(n.n2))
    N2′ = kronecker(𝕀(2), 𝕀(n.n1), N(n.n2))
    for i in 1:length(ϵ)
        if (ϵ[i] - E_ground) < 1e-16
            n1 += tr(N1′*P[:,i]*adjoint(P[:,i]))
            n2 += tr(N2′*P[:,i]*adjoint(P[:,i]))
        end
    end
    return (n1 + n2)/(n.n1 + n.n2)
end
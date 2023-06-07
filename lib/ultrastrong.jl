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

function nTzero(prm::Lorentzian, ang::CouplingAngles, n::Levels) # Ground-state number expectation
    H = HTot(prm, ang, n)
    F = eigen(H)
    ϵ, P = F.values, F.vectors
    E_ground = ϵ[1]
    n1 = 0
    n2 = 0
    N1 = kronecker(𝕀(2), N(n.n1), 𝕀(n.n2))
    N2 = kronecker(𝕀(2), 𝕀(n.n1), N(n.n2))
    for i in 1:length(ϵ)
        if (ϵ[i] - E_ground) < 1e-16
            n1 += tr(N1*P[:,i]*adjoint(P[:,i]))
            n2 += tr(N2*P[:,i]*adjoint(P[:,i]))
        end
    end
    return (n1 + n2)/2
end

function nTzero′(prm::Lorentzian, ang::CouplingAngles, n::Levels) # Ground-state number expectation
    H = HTot′(prm, ang, n)
    F = eigen(H)
    ϵ, P = F.values, F.vectors
    E_ground = ϵ[1]
    n1 = 0
    n2 = 0
    N1 = kronecker(𝕀(2), N(n.n1), 𝕀(n.n2))
    N2 = kronecker(𝕀(2), 𝕀(n.n1), N(n.n2))
    for i in 1:length(ϵ)
        if (ϵ[i] - E_ground) < 1e-16
            n1 += tr(N1*P[:,i]*adjoint(P[:,i]))
            n2 += tr(N2*P[:,i]*adjoint(P[:,i]))
        end
    end
    return (n1 + n2)/2
end

function osc_positions(prm, ang, n)
    probabilities = diag(ρGroundBath(prm, ang, n))
    X1(n) = kronecker(𝕀(2), X(n.n1), 𝕀(n.n2))
    X2(n) = kronecker(𝕀(2), 𝕀(n.n1), X(n.n2))
    xmax = 1/sqrt(2)*maximum(eigvals(X(n.n1)))
    xvals = LinRange(-1.5*xmax, 1.5*xmax, 500)
    nvals = 1:length(probabilities)
    dist = zeros(length(nvals), length(xvals))
    for i in 1:length(nvals)
        basis = FockBasis(n.n1)
        psi = fockstate(basis, nvals[i])
        dist[i,:] = probabilities[i]*sum(wigner(psi, xvals, xvals), dims=1)
    end
    return xvals, sum(dist, dims=1)
end
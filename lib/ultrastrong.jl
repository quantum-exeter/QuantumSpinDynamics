########################
#### ultrastrong.jl ####
########################

#### T=0 sz magnetisation for UNTRANSFORMED Hamiltonian ####
function szTzero(prm::Lorentzian, ang::CouplingAngles, n::Levels)
    H = HTot(prm, ang, n) 
    F = eigen(H)
    ϵ, P = F.values, F.vectors
    E_ground = ϵ[1]
    sz = 0
    Sz = kronecker(σz, 𝕀(n.n1), 𝕀(n.n2))
    for i in 1:length(ϵ)
        if (ϵ[i] - E_ground) < 1e-16
            sz += tr(Sz*P[:,i]*adjoint(P[:,i]))
        end
    end
    return sz
end

#### T=0 sz magnetisation for TRANSFORMED Hamiltonian ####
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

#### T=0 average number expectation for UNTRANSFORMED Hamiltonian ####
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

#### T=0 average number expectation for TRANSFORMED Hamiltonian ####
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
    return n1, n2
end

xmax(n) = 1/sqrt(2)*maximum(eigvals(X(n.n1)))

function osc_positions(prm, ang, n)
    probabilities = diag(ρGroundBath(prm, ang, n)) # Generates the probabilities for occupation of the oscillator levels
    xvals = collect(LinRange(-1.5*xmax(n), 1.5*xmax(n), 500))  # Position sweep values
    nvals = 1:length(probabilities) # All the possible oscillator levels
    dist = zeros(length(xvals), length(nvals)) # Set up empty array for position distribution
    for i in 1:length(nvals)
        basis = FockBasis(n.n1) 
        psi = fockstate(basis, nvals[i])
        dist[:,i] = probabilities[i]*sum(wigner(psi, xvals, xvals), dims=1)
    end
    return hcat(xvals, sum(dist, dims=2))
end
###################
#### states.jl ####
###################

### Thermal State ###
function thermal(H, T)
    n = size(H, 1)
    ϵ = eigen(H).values
    P = eigen(H).vectors
    𝒵 = sum(exp(-ϵ[i]/T) for i = 1:n)
    ρ = (1/𝒵)*Diagonal([exp(-ϵ[i]/T) for i = 1:n])
    return P*ρ*adjoint(P)
end

ρGibbs(T) = thermal(HS(), T)

### MFGS ###
function ρMFGS(prm::Lorentzian, ang::CouplingAngles, n::Levels, T)
    H = HTot(prm, ang, n) # Hamiltonian for the given coupling dimension
    proj = eigen(H).vectors # Projector onto the basis of H
    HTr = adjoint(proj)*H*proj # Transformed H
    stateTot = proj*thermal(HTr, T)*adjoint(proj)
    ntr = Int(hspace_size(n)/2)
    return ptrace(stateTot, ntr)
end

### Ground States ###
function ρGround(prm::Lorentzian, ang::CouplingAngles, n::Levels)
    H = HTot(prm, ang, n)
    state = eigen(H).vectors[:,1]
    return adjoint(state)*state
end

### Initial States (for Dynamics) ###

## Bloch State ##
function ρBloch()
    α = -π/2
    β = 0
    return [cos(α/2)^2 0.5*exp(-im*β)*sin(α); 0.5*exp(im*β)*sin(α) sin(α/2)^2]
end

## Joint States ##
ρ0(prm::LorPrm1D, n::Lev1D, T) = kronecker(ρBloch(), thermal(prm.ω01*N(n.n1), T))
ρ0(prm::LorPrm2D, n::Lev2D, T) = kronecker(ρBloch(), thermal(prm.ω01*N(n.n1), T), thermal(prm.ω02*N(n.n2), T))
ρ0(prm::LorPrm3D, n::Lev3D, T) = kronecker(ρBloch(), thermal(prm.ω01*N(n.n1), T), thermal(prm.ω02*N(n.n2), T), thermal(prm.ω03*N(n.n3), T))
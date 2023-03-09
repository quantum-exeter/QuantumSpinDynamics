###################
#### states.jl ####
###################

### Thermal State ###
function thermal(H, T)
    n = size(H, 1)
    F = eigen(H)
    ϵ, P = F.values, F.vectors
    # ϵ = eigen(H).values
    # P = eigen(H).vectors
    𝒵 = sum(exp(-ϵ[i]/T) for i = 1:n)
    ρ = (1/𝒵)*Diagonal([exp(-ϵ[i]/T) for i = 1:n])
    return P*ρ*adjoint(P)
end

ρGibbs(T) = thermal(HS(), T)

### MFGS ###
function ρMFGS(prm::Lorentzian, ang::CouplingAngles, n::Levels, T)
    H = HTot(prm, ang, n) # Hamiltonian for the given coupling dimension
    # proj = eigen(H).vectors # Projector onto the basis of H
    # HTr = adjoint(proj)*H*proj # Transformed H
    # stateTot = proj*thermal(HTr, T)*adjoint(proj)
    stateTot = thermal(H, T)
    ntr = Int(hspace_size(n)/2)
    rho = ptrace(stateTot, ntr)
    return rho
end

### Variational States ###
function ψGround(prm::LorPrm3D, ang::CouplAng3D, n::Lev3D)
    H = HTot(prm, ang, n)
    state = eigen(H).vectors[:,1]
    return state
end

function ϕp(prm::LorPrm3D, ang::CouplAng3D, n::Lev3D)
    su = [1 0]
    state = chopBoth.(kronecker(su, 𝕀(n.n1*n.n2*n.n3))*ψGround(prm, ang, n))
    stateNorm = (1/norm(state))*state
    return stateNorm
end

function ϕm(prm::LorPrm3D, ang::CouplAng3D, n::Lev3D)
    sd = [0 1]
    state = chopBoth.(kronecker(sd, 𝕀(n.n1*n.n2*n.n3))*ψGround(prm, ang, n))
    stateNorm = (1/norm(state))*state
    return stateNorm
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
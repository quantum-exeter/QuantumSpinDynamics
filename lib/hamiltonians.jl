#########################
#### hamiltonians.jl ####
#########################

### Creation and annihilation operators ###

## Creation operator ##
function create(n)
    matrix = zeros(n, n)
    for i in 1:n
        new_row = zeros(1, n)
        for j in Array(1:n)
            if i == j+1
                new_row[j] = sqrt(i-1)
            else
                new_row[j] = 0
            end
        end
        matrix[[i],:] = new_row
    end
    return matrix
end

## Annihilation operator ##
annihilate(n) = adjoint(create(n))

## Position operator ##
X(n) = annihilate(n) + create(n)
N(n) = create(n)*annihilate(n)

### Hamiltonians ###

## Spin Hamiltonian ##
HS() = -sign(γ)*σz
HS(n::Lev1D) = -sign(γ)*kronecker(σz, 𝕀(n.n1))
HS(n::Lev2D) = -sign(γ)*kronecker(σz, 𝕀(n.n1), 𝕀(n.n2))
HS(n::Lev3D) = -sign(γ)*kronecker(σz, 𝕀(n.n1), 𝕀(n.n2), 𝕀(n.n3))

## Interaction Hamiltonian ##
HInt(prm::LorPrm1D, ang::CouplAng1D, n::Lev1D) = sqrt(3)*sqrt(prm.α1/prm.ω01)*kronecker(sc(ang.θ1, ang.ϕ1), X(n.n1))
HInt(prm::LorPrm2D, ang::CouplAng2D, n::Lev2D) = sqrt(prm.α1/prm.ω01)*kronecker(sc(ang.θ1, ang.ϕ1), X(n.n1), 𝕀(n.n2)) + sqrt(prm.α2/prm.ω02)*kronecker(sc(ang.θ2, ang.ϕ2), 𝕀(n.n1), X(n.n2))
HInt(prm::LorPrm3D, ang::CouplAng3D, n::Lev3D) = sqrt(prm.α1/prm.ω01)*kronecker(sc(ang.θ1, ang.ϕ1), X(n.n1), 𝕀(n.n2), 𝕀(n.n3)) + sqrt(prm.α2/prm.ω02)*kronecker(sc(ang.θ2, ang.ϕ2), 𝕀(n.n1), X(n.n2), 𝕀(n.n3)) + sqrt(prm.α3/prm.ω03)*kronecker(sc(ang.θ3, ang.ϕ3), 𝕀(n.n1), 𝕀(n.n2), X(n.n3))

## Bath Hamiltonian ##
HB(prm::LorPrm1D, n::Lev1D) = (1/s0)*(prm.ω01*kronecker(𝕀s, N(n.n1)))
HB(prm::LorPrm2D, n::Lev2D) = (1/s0)*(prm.ω01*kronecker(𝕀s, N(n.n1), 𝕀(n.n2)) + prm.ω02*kronecker(𝕀s, 𝕀(n.n1), N(n.n2)))
HB(prm::LorPrm3D, n::Lev3D) = (1/s0)*(prm.ω01*kronecker(𝕀s, N(n.n1), 𝕀(n.n2), 𝕀(n.n3)) + prm.ω02*kronecker(𝕀s, 𝕀(n.n1), N(n.n2), 𝕀(n.n3)) + prm.ω03*kronecker(𝕀s, 𝕀(n.n1), 𝕀(n.n2), N(n.n3)))

## Total Hamiltonian ##
HTot(prm::Lorentzian, ang::CouplingAngles, n::Levels) = HS(n) + HInt(prm, ang, n) + HB(prm, n) 
#########################
#### hamiltonians.jl ####
#########################

### Creation and Annihilation Operators ###

## Creation Operator ##
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

## Annihilation Operator ##
function annihilate(n)
    return adjoint(create(n))
end

### Hamiltonians ###

## Bare Spin Hamiltonian ##
HG() = -sign(γ)*σz

## 1D Hamiltonian ##
function HS1D(n1, λ1, Ω1)
    spin = -sign(γ)*kronecker(σz, 𝕀(n1))
    rc = scale*kronecker(𝕀s, (Ω1/ωL)*(create(n1)*annihilate(n1)))
    int = (λ1/ωL)*kronecker(sc(1), (create(n1) + annihilate(n1)))
    return spin + rc + int
end

## 2D Hamiltonian ##
function HS2D(n1, n2, λ1, λ2, Ω1, Ω2)
    spin = -sign(γ)*kronecker(σz, 𝕀(n1), 𝕀(n2))
    rc = scale*kronecker(𝕀s, (Ω1/ωL)*(create(n1)*annihilate(n1)), 𝕀(n2)) + scale*kronecker(𝕀s, 𝕀(n1), (Ω2/ωL)*(create(n2)*annihilate(n2)))
    int = (λ1/ωL)*kronecker(sc(1), (create(n1) + annihilate(n1)), 𝕀(n2)) + (λ2/ωL)*kronecker(sc(2), 𝕀(n1), (create(n2) + annihilate(n2)))
    return spin + rc + int
end

## 3D Hamiltonian ##
function HS3D(n1, n2, n3, λ1, λ2, λ3, Ω1, Ω2, Ω3)
    spin = -sign(γ)*kronecker(σz, 𝕀(n1), 𝕀(n2), 𝕀(n3))
    rc = scale*kronecker(𝕀s, (Ω1/ωL)*(create(n1)*annihilate(n1)), 𝕀(n2), 𝕀(n3)) + scale*kronecker(𝕀s, 𝕀(n1), (Ω2/ωL)*(create(n2)*annihilate(n2)), 𝕀(n3)) + scale*kronecker(𝕀s, 𝕀(n1), 𝕀(n2), (Ω3/ωL)*(create(n3)*annihilate(n3)))
    int = (λ1/ωL)*kronecker(sc(1), (create(n1) + annihilate(n1)), 𝕀(n2), 𝕀(n3)) + (λ2/ωL)*kronecker(sc(2), 𝕀(n1), (create(n2) + annihilate(n2)), 𝕀(n3)) + (λ3/ωL)*kronecker(sc(3), 𝕀(n1), 𝕀(n2), (create(n3) + annihilate(n3)))
    return spin + rc + int
end

## Single Function for Access to All Dimensions ##
function HS(dim)
    if dim == 1
        return HS1D(n1, λ1, Ω1)
    elseif dim == 2
        return HS2D(n1, n2, λ1, λ2, Ω1, Ω2)
    elseif dim == 3 
        return HS3D(n1, n2, n3, λ1, λ2, λ3, Ω1, Ω2, Ω3)
    else
        print("Please return a dimension of either 1, 2 or 3.")
    end
end

## General Bath Hamiltonian ##
HB(n, Ω) = (Ω/ωL)*(create(n)*annihilate(n))
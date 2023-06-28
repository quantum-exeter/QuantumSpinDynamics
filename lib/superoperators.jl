###########################
#### superoperators.jl ####
###########################

### Superoperator for the Redfield equation ###
function 𝒮(prm::Lorentzian, ang::CouplingAngles, n::Levels, T)

    H = HTot(prm, ang, n) # system Hamiltonian for given dimension of coupling

    ℒ(operator) = kronecker(operator, 𝕀(hspace_size(n))) # define the left multiplication superoperator
    ℛ(operator) = kronecker(𝕀(hspace_size(n)), transpose(operator)) # define the right multiplication superoperator

    supop = -im*(ℒ(H) - ℛ(H))

    for i in 1:dim(n)
        trans = transitions(prm, ang, n, i)
        ωb = trans[1]
        ATr = trans[2]
        len = length(ωb)
        ATot = sum(ATr[j] for j = 1:len)
        χ = (π/2)*sum(spectral_density_Ohm(ωb[j], prm)[i]*coth((ωb[j])/(2*T))*ATr[j] for j = 1:len) 
        Θ = (π/2)*sum(spectral_density_Ohm(ωb[j], prm)[i]*ATr[j] for j = 1:len)
        supop += - ℒ(ATot)*(ℒ(χ) - ℛ(χ)) + ℛ(ATot)*(ℒ(χ) - ℛ(χ)) + ℒ(ATot)*(ℒ(Θ) + ℛ(Θ)) - ℛ(ATot)*(ℒ(Θ) + ℛ(Θ))
    end

    return supop

end
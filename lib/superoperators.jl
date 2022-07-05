###########################
#### superoperators.jl ####
###########################

### Superoperator for the Redfield Equation ###
function 𝒮(prm::Lorentzian, ang::CouplingAngles, n::Levels, T)

    ## Useful Definitions ##
    H = HTot(prm, ang, n) # System Hamiltonian for given dimension of coupling

    ## Bohr Frequencies/Jump Operators for Each Coupling Dimension ##
    ωb(i) = transitions(prm, ang, n, i)[1] # Rewrite function outputs (Bohr freqs) in more compact form
    len(i) = length(ωb(i)) # Find the length of list (for each bath) to iterate over
    ATr(i) = transitions(prm, ang, n, i)[2] # Rewrite function outputs (transformed jump ops) in more compact form
    ATot(i) = sum(ATr(i)[j] for j = 1:len(i)) # Define the sum of all jump operators for each bath

    ## Iles-Smith Superoperators ##
    χ(i) = (π/2)*sum(spectral_density(ωb(i)[j], prm)[i]*coth((ωb(i)[j])/(2*T))*ATr(i)[j] for j = 1:len(i)) # Iles-Smith χ superoperator
    Θ(i) = (π/2)*sum(spectral_density(ωb(i)[j], prm)[i]*ATr(i)[j] for j = 1:len(i)) # Iles-Smith Θ superoperator

    ## Left/Right Multiplication Superoperators ##
    ℒ(operator) = kronecker(operator, 𝕀(hspace_size(n))) # Define the left multiplication superoperator
    ℛ(operator) = kronecker(𝕀(hspace_size(n)), transpose(operator)) # Define the right multiplication superoperator

    ## Return the Superoperator ##
    return -im*(ℒ(H) - ℛ(H)) + sum(- ℒ(ATot(i))*(ℒ(χ(i)) - ℛ(χ(i))) + ℛ(ATot(i))*(ℒ(χ(i)) - ℛ(χ(i))) + ℒ(ATot(i))*(ℒ(Θ(i)) + ℛ(Θ(i))) - ℛ(ATot(i))*(ℒ(Θ(i)) + ℛ(Θ(i))) for i in 1:dim(n))

end
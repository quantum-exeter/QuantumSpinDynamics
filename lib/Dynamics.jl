module Dynamics

    ### Import Packages ###
    using LinearAlgebra
    using Kronecker

    ### Inclusions ###
    include("constants.jl")
    include("variables.jl")
    include("maths.jl")
    include("initialStates.jl")
    include("hamiltonians.jl")
    include("spectralDensity.jl")

    ####################################
    ####################################
    ####################################

    ### Exports ###
    export 𝒮, ρ0, hspace_dimension, σx, σy, σz, 𝕀, realIfClose
    
    ### Dynamics-Specific Functions ###

    ## Transition Frequencies and Jump Operators ##

    # Dimensionless Position Operator (Jump Function for this System) #
    X(n) = create(n) + annihilate(n)

    # 1D Coupling #
    jump1D(n1) = [kronecker(𝕀s, X(n1))]

    # 2D Coupling #
    jump2D(n1, n2) = [kronecker(𝕀s, X(n1), 𝕀(n2)), kronecker(𝕀s, 𝕀(n1), X(n2))]

    # 3D Coupling #
    jump3D(n1, n2, n3) = [kronecker(𝕀s, X(n1), 𝕀(n2), 𝕀(n3)), kronecker(𝕀s, 𝕀(n1), X(n2), 𝕀(n3)), kronecker(𝕀s, 𝕀(n1), 𝕀(n2), X(n3))] 

    # Single Function for Access to All Dimensions #
    function jump(dim)
        if dim == 1
            return jump1D(n1)
        elseif dim == 2
            return jump2D(n1, n2)
        elseif dim == 3
            return jump3D(n1, n2, n3)
        else
            print("Please return a dimension of either 1, 2 or 3.")
        end
    end

    function transitions(dim, i) # dim is the coupling dimension (1 for 1D, 2 for 2D etc...) and i the specific RC to generate Bohr frequencies/jump operators for
        H = HS(dim)
        A = jump(dim)[i]
        n = hspace_dimension(dim)
        table = zeros(n, n)
        bohr_freqs = Float64[]
        jump_ops = Any[]
        eval = eigen(H).values
        evec = eigen(H).vectors
        proj(j) = evec[:,j]*adjoint(evec[:,j])
        # Create a table of transition frequencies
        for k in 1:n
            for l in 1:n
                table[k, l] = eval[l] - eval[k]
                if table[k, l] !=  0.0
                    append!(bohr_freqs, table[k, l])
                    push!(jump_ops, proj(l)*A*proj(k))
                end
            end
        end
        return bohr_freqs, jump_ops
    end

    ## Superoperator for the Redfield Equation ##
    function 𝒮(dim)

        ## Useful Definitions ##
        H = HS(dim) # System Hamiltonian for given dimension of coupling
        A = jump(dim) # Jump operator list for given dimension of coupling in ORIGINAL basis
        n = hspace_dimension(dim) # Dimension of Hilbert space for given coupling dimension

        ## Bohr Frequencies/Jump Operators for Each Coupling Dimension ##
        transitions_list(i) = transitions(dim, i) # Extract Bohr freq-jump op pairs for the ith bath
        len(i) = length(transitions_list(i)[1]) # Find the length of list (for each bath) to iterate over
        ωb(i) = transitions_list(i)[1] # Rewrite function outputs (Bohr freqs) in more compact form
        ATr(i) = transitions_list(i)[2] # Rewrite function outputs (transformed jump ops) in more compact form
        ATot(i) = sum(ATr(i)[j] for j = 1:len(i)) # Define the sum of all jump operators for each bath

        ## Iles-Smith Superoperators ##
        χ(i) = (π/2)*sum(spectral_density(ωb(i)[j], δ_list(i))*coth((ωb(i)[j])/(2*TDyn))*ATr(i)[j] for j = 1:len(i)) # Iles-Smith χ superoperator
        Θ(i) = (π/2)*sum(spectral_density(ωb(i)[j], δ_list(i))*ATr(i)[j] for j = 1:len(i)) # Iles-Smith Θ superoperator

        ## Left/Right Multiplication Superoperators ##
        ℒ(operator) = kronecker(operator, 𝕀(n)) # Define the left multiplication superoperator
        ℛ(operator) = kronecker(𝕀(n), transpose(operator)) # Define the right multiplication superoperator

        ## Return the Superoperator ##
        return sum(-im*(ℒ(H) - ℛ(H)) - ℒ(ATot(i))*(ℒ(χ(i)) - ℛ(χ(i))) + ℛ(ATot(i))*(ℒ(χ(i)) - ℛ(χ(i))) + ℒ(ATot(i))*(ℒ(Θ(i)) + ℛ(Θ(i))) - ℛ(ATot(i))*(ℒ(Θ(i)) + ℛ(Θ(i))) for i in 1:dim)

    end

end
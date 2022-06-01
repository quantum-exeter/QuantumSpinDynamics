module Dynamics

    ### Import Modules ###
    using LinearAlgebra
    using Kronecker
    using OrdinaryDiffEq
    using DifferentialEquations

    ####################################
    ####################################
    ####################################

    ### Exports ###
    export 𝒮

    ### Dynamics-Specific Functions ###

    ## Transition Frequencies ##
    function transitions(n, λ, Ω)
        table = zeros(2n, 2n)
        bohr_frequencies = Float64[]
        jump_ops = Any[]
        eval = eigvals(HRC1D(n, λ, Ω))
        evec = eigvecs(HRC1D(n, λ, Ω))
        proj(i) = evec[:,i]*adjoint(evec[:,i])
        A = kron(𝕀s, (create(n) + annihilate(n)))
        # Create a table of transition frequencies
        for i in 1:2n
            for j in 1:2n
                table[i,j] = eval[j] - eval[i]
                if table[i,j] !=  0.0
                    append!(bohr_frequencies, table[i,j])
                    push!(jump_ops, proj(j)*A*proj(i))
                end
            end
        end
        return(bohr_frequencies, jump_ops)
    end

    ## Iles-Smith RC ME Parameters ##
    function χop(n, Ω, λ, δ, T)
        transitions_list = transitions(n, λ, Ω)
        len = length(transitions_list[1])
        ωb = transitions_list[1]
        Aj = transitions_list[2]
        return((pi/2)*sum(spectral_density(ωb[i], δ)*coth((cfac*ωb[i])/(2*T))*Aj[i] for i=1:len))
    end

    function Θop(n, Ω, λ, δ)
        transitions_list = transitions(n, λ, Ω)
        len = length(transitions_list[1])
        ωb = transitions_list[1]
        Aj = transitions_list[2]
        return((π/2)*sum(spectral_density(ωb[i], δ)*Aj[i] for i=1:len))
    end

    ## Left/Right Multiplication Superoperators ##
    𝕀(n) = 1*Matrix(I, 2n, 2n)
    ℒ(operator, n) = kronecker(operator, 𝕀(n))
    ℛ(operator, n) = kronecker(𝕀(n), transpose(operator))

    function 𝒮(n, Ω, λ, δ, T)
        L(operator) = ℒ(operator, n)
        R(operator) = ℛ(operator, n)
        H = HRC1D(n, λ, Ω)
        transitions_list = transitions(n, λ, Ω)
        len = length(transitions_list[1])
        Atot = sum(transitions_list[2][i] for i=1:len)
        χ = χop(n, Ω, λ, δ, T)
        Θ = Θop(n, Ω, λ, δ)
        return(-im*(L(H) - R(H)) - L(Atot)*(L(χ) - R(χ)) + R(Atot)*(L(χ) - R(χ)) + L(Atot)*(L(Θ) + R(Θ)) - R(Atot)*(L(Θ) + R(Θ)))
    end
    
end
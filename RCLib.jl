module RCLib

    ### Import Modules ###
    using LinearAlgebra
    using Kronecker
    using OrdinaryDiffEq
    using DifferentialEquations
    using Plots

    ####################################
    ####################################
    ####################################

    ### Exports ###
    export ωL, ρ0, 𝒮, sx0, sy0, sz0, 𝕀b

    ### Variables ###
    γ = -1.76*10^(11) # Gyromagnetic ratio for an electron (T^-1s^-1)
    Bext = 10 # External magnetic field (T)
    ωL = abs(γ)*Bext # Larmor frequency (s^-1)
    ħ = 1.05*10^(-34) # Reduced Planck constant (Js)
    kB = 1.38*10^(-23) # Boltzmann cosntant (JK^-1)
    Λ = 10^10 # Cutoff frequency for spectral density
    cfac = (ħ*ωL)/kB

    # Identity matrices for subsystem and bath #
    𝕀s = 1*Matrix(I, 2, 2)
    𝕀b(n) = 1*Matrix(I, n, n)

    ### General Maths Functions ###

    # Commutators/Anticommutators #
    comm(A,B) = A*B - B*A
    acomm(A,B) = A*B + B*A

    ### Initial States ###

    # Spins #
    σx = [[0 1];[1 0]]
    σy = [[0 -im];[im 0]]
    σz = [[1 0];[0 -1]]
    sx0 = 0.5*σx
    sy0 = 0.5*σy
    sz0 = 0.5*σz

    function bloch_state(θ, ϕ)
        return([cos(θ/2)^2 0.5*exp(-im*ϕ)*sin(θ); 0.5*exp(im*ϕ)*sin(θ) sin(θ/2)^2])
    end

    # RC #
    function thermal_state(n, ω0, T)
        ϵ = zeros(n) # Initialise empty array for energies of RC harmonic oscillator
        for i in 1:n
            ϵ[i] = ħ*ω0*(i - 0.5)
        end
        𝒵 = sum(exp(-ϵ[i]/(kB*T)) for i = 1:n)
        Hb = diagm(ϵ)
        return((1/𝒵)*exp(-Hb/(kB*T)))
    end

    # Joint Initital State #
    ρ0(θ, ϕ, n, ω0, T) = kronecker(bloch_state(θ, ϕ), thermal_state(n, ω0, T))

    ### Creation and Annihilation Operators ###
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

    function annihilate(n)
        matrix = zeros(n, n)
        for i in Array(1:n)
            new_row = zeros(1, n)
            for j in 1:n
                if i == j+1
                    new_row[j] = sqrt(i-1)
                else
                    new_row[j] = 0
                end
            end
            matrix[[i],:] = new_row
        end
        return adjoint(matrix)
    end

    ### RC Hamiltonian ###
    HRC(n, λ, Ω) = kronecker(sz0, 𝕀b(n)) + (λ/ωL)*kronecker(sx0, (create(n) + annihilate(n))) + kronecker(𝕀s, (Ω/ωL)*(create(n)*annihilate(n)))

    ### Transition Frequencies ###
    function transitions(n, λ, Ω)
        table = zeros(2n, 2n)
        bohr_frequencies = Float64[]
        jump_ops = Any[]
        eval = eigvals(HRC(n, λ, Ω))
        evec = eigvecs(HRC(n, λ, Ω))
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

    ### Spectral Density ###
    spectral_density(ω, δ) = δ*ω*(Λ^2/(ω^2 + Λ^2))

    ### Redfield (RF) Equation ###

    # Iles-Smith RC ME Parameters
    function χop(n, Ω, λ, δ, T)
        len = length(transitions(n, λ, Ω)[1])
        ωb = transitions(n, λ, Ω)[1]
        Aj = transitions(n, λ, Ω)[2]
        return((pi/2)*sum(spectral_density(ωb[i], δ)*coth((cfac*ωb[i])/(2*T))*Aj[i] for i=1:len))
    end

    function Θop(n, Ω, λ, δ)
        len = length(transitions(n, λ, Ω)[1])
        ωb = transitions(n, λ, Ω)[1]
        Aj = transitions(n, λ, Ω)[2]
        return(-(π/2)*sum(spectral_density(ωb[i], δ)*Aj[i] for i=1:len))
    end

    # Left/Right Multiplication Superoperators
    𝕀(n) = 1*Matrix(I, 2n, 2n)
    ℒ(operator, n) = kronecker(operator, 𝕀(n))
    ℛ(operator, n) = kronecker(𝕀(n), transpose(operator))

    function 𝒮(n, Ω, λ, δ, T)
        L(operator) = ℒ(operator, n)
        R(operator) = ℛ(operator, n)
        H = HRC(n, λ, Ω)
        len = length(transitions(n, λ, Ω)[1])
        Atot = sum(transitions(n, λ, Ω)[2][i] for i=1:len)
        χ = χop(n, Ω, λ, δ, T)
        Θ = Θop(n, Ω, λ, δ)
        return(-im*(L(H) - R(H)) - L(Atot)*(L(χ) - R(χ)) + R(Atot)*(L(χ) - R(χ)) + L(Atot)*(L(Θ) + R(Θ)) - R(Atot)*(L(Θ) + R(Θ)))
    end
end

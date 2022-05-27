module RCLib

    ### Import Modules ###
    using LinearAlgebra
    using Kronecker
    using OrdinaryDiffEq
    using DifferentialEquations

    ####################################
    ####################################
    ####################################

    ### Exports ###
    export ωL, real_if_close, ρ0, 𝒮, sx0, sy0, sz0, 𝕀b, gibbs, HSpG, HRC1D, HRC2D, HRC3D, ptrace, ℱ, pred

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

    # Square
    square(n) = n * n

    # Make Integer
    int(x) = floor(Int, x)

    # Check Size Im Parts
    real_if_close(c) = isnan(imag(c)) || imag(c) < 1e-14 ? real(c) : c
    real_if_close(c::AbstractArray) = real_if_close.(c)

    ### Initial States ###

    # Spins #
    σx = [[0 1];[1 0]]
    σy = [[0 -im];[im 0]]
    σz = [[1 0];[0 -1]]
    sx0 = 0.5*σx
    sy0 = 0.5*σy
    sz0 = 0.5*σz

    function bloch_state(θ, ϕ)
        return [cos(θ/2)^2 0.5*exp(-im*ϕ)*sin(θ); 0.5*exp(im*ϕ)*sin(θ) sin(θ/2)^2]
    end

    # RC #
    function gibbs(H, T)
        n = size(H,1)
        ϵ = eigen(H).values
        P = eigen(H).vectors
        𝒵 = sum(exp(-(cfac*ϵ[i])/T) for i = 1:n)
        ρ = (1/𝒵)*Diagonal([exp(-(cfac*ϵ[i])/T) for i = 1:n])
        return P*ρ*adjoint(P)
    end

    # Joint Initital State #
    function ρ0(θ, ϕ, Ω, n, T)
        H_bath = ((Ω/ωL)*(create(n)*annihilate(n)))
        thermal_state = gibbs(H_bath, T)
        return kronecker(bloch_state(θ, ϕ), thermal_state)
    end

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
        return adjoint(create(n))
    end

    ### Gibbs Spin Hamiltonian ###
    HSpG = -sign(γ)*sz0

    ### 1D RC Hamiltonian ###
    HRC1D(n, λ, Ω) = -sign(γ)*kronecker(sz0, 𝕀b(n)) + (λ/ωL)*kronecker(sx0, (create(n) + annihilate(n))) + kronecker(𝕀s, (Ω/ωL)*(create(n)*annihilate(n)))

    ### 2D RC Hamiltonian ###
    function HRC2D(nx, nz, λx, λz, Ωx, Ωz)
        spin = -sign(γ)*kronecker(sz0, 𝕀b(nx), 𝕀b(nz))
        rc = kronecker(𝕀s, (Ωx/ωL)*(create(nx)*annihilate(nx)), 𝕀b(nz)) + kronecker(𝕀s, 𝕀b(nx), (Ωz/ωL)*(create(nz)*annihilate(nz)))
        int = (λx/ωL)*kronecker(sx0, (create(nx) + annihilate(nx)), 𝕀b(nz)) + (λz/ωL)*kronecker(sz0, 𝕀b(nx), (create(nz) + annihilate(nz)))
        return(spin + rc + int)
    end

    ### 3D RC Hamiltonian ###
    function HRC3D(nx, ny, nz, λx, λy, λz, Ωx, Ωy, Ωz)
        spin = -sign(γ)*kronecker(sz0, 𝕀b(nx), 𝕀b(ny), 𝕀b(nz))
        rc = kronecker(𝕀s, (Ωx/ωL)*(create(nx)*annihilate(nx)), 𝕀b(ny), 𝕀b(nz)) + kronecker(𝕀s, 𝕀b(nx), (Ωy/ωL)*(create(ny)*annihilate(ny)), 𝕀b(nz)) + kronecker(𝕀s, 𝕀b(nx), 𝕀b(ny), (Ωz/ωL)*(create(nz)*annihilate(nz)))
        int = (λx/ωL)*kronecker(sx0, (create(nx) + annihilate(nx)), 𝕀b(ny), 𝕀b(nz)) + (λy/ωL)*kronecker(sy0, 𝕀b(nx), (create(ny) + annihilate(ny)), 𝕀b(nz)) + (λz/ωL)*kronecker(sz0, 𝕀b(nx), 𝕀b(ny), (create(nz) + annihilate(nz)))
        return(spin + rc + int)
    end

    ### Transition Frequencies ###
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

    ### Spectral Density ###
    spectral_density(ω, δ) = δ*ω*(Λ^2/(ω^2 + Λ^2))

    ### Dynamics ###

    # Iles-Smith RC ME Parameters
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

    # Left/Right Multiplication Superoperators
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

    ### HMF Calculations ###

    # Partial Trace
    function ptrace(ρ, n)
        nR = int(size(ρ, 1)/n)
        return(sum((𝕀b(nR)⊗(𝕀b(n)[[i],:]))*ρ*(𝕀b(nR)⊗(𝕀b(n)[:,i])) for i=1:n))
    end

    # Uhlmann Fidelity
    ℱ(ρ1, ρ2) = square(tr(sqrt(sqrt(ρ1)*ρ2*sqrt(ρ1))))

    # Classical 3D Gibbs/MFGS S_z Function
    pred(T) = coth(cfac/(2*T)) - (2*T)/cfac

end
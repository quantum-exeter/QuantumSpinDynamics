########################
#### transitions.jl ####
########################

### Jump Operators ###
jump(n::Lev1D) = [kronecker(𝕀s, X(n.n1))]
jump(n::Lev2D) = [kronecker(𝕀s, X(n.n1), 𝕀(n.n2)), kronecker(𝕀s, 𝕀(n.n1), X(n.n2))]
jump(n::Lev3D) = [kronecker(𝕀s, X(n.n1), 𝕀(n.n2), 𝕀(n.n3)), kronecker(𝕀s, 𝕀(n.n1), X(n.n2), 𝕀(n.n3)), kronecker(𝕀s, 𝕀(n.n1), 𝕀(n.n2), X(n.n3))] 

### Function that generates Bohr freq/jump operator pairings for the ith coupling ###
function transitions(prm::Lorentzian, ang::CouplingAngles, n::Levels, i)
    H = HTot(prm, ang, n)
    A = jump(n)[i]
    d = hspace_size(n)
    table = zeros(d, d)
    bohr_freqs = Float64[]
    jump_ops = Any[]
    eval = eigen(H).values
    evec = eigen(H).vectors
    proj(j) = evec[:,j]*adjoint(evec[:,j])
    # Create a table of transition frequencies
    for k in 1:d
        for l in 1:d
            table[k, l] = eval[l] - eval[k]
            if table[k, l] !=  0.0
                append!(bohr_freqs, table[k, l])
                push!(jump_ops, proj(l)*A*proj(k))
            end
        end
    end
    return bohr_freqs, jump_ops
end
alg = Vern7()

function dsolve(prm::Lorentzian, ang::CouplingAngles, n::Levels, T, tspan)
    superop = 𝒮(prm, ang, n, T)
    state_init = ρ0(prm, n, T)
    dstate(dρ, ρ, v, t) = mul!(dρ, superop, ρ) # Solves the DE
    prob = ODEProblem(dstate, vec(state_init), tspan)
    solution = solve(prob, alg)
    out(t) = reshape(solution(t), (hspace_size(n), hspace_size(n))) # Reformats vector into a density matrix
    return out
end

#########################
#### diffEqSolver.jl ####
#########################

### Algorithm ###
alg = Vern7()

### Differential equation solver ###
function dsolve(prm::Lorentzian, ang::CouplingAngles, n::Levels, T, tspan)
    superop = 𝒮(prm, ang, n, T)
    state_init = ρ0(prm, n, T) # initial state
    dstate(dρ, ρ, v, t) = mul!(dρ, superop, ρ) # solves the DE
    prob = ODEProblem(dstate, vec(state_init), tspan)
    solution = solve(prob, alg)
    out(t) = reshape(solution(t), (hspace_size(n), hspace_size(n))) # reformats vector into a density matrix
    return out
end

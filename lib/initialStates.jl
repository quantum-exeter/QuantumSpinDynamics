# Spins #
σx = [[0 1];[1 0]]
σy = [[0 -im];[im 0]]
σz = [[1 0];[0 -1]]
sx0 = 0.5*σx
sy0 = 0.5*σy
sz0 = 0.5*σz

function s0(i)
    θ = [π/2 0 π/2]
    ϕ = [0 0 π/2]
    return sx0*(sin(θ[i])*cos(ϕ[i])) + sy0*(sin(θ[i])*sin(ϕ[i])) + sz0*cos(θ[i])
end

function bloch_state(α, β)
    return [cos(α/2)^2 0.5*exp(-im*β)*sin(α); 0.5*exp(im*β)*sin(α) sin(α/2)^2]
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
function ρ0(α, β, Ω, n, T)
    H_bath = ((Ω/ωL)*(create(n)*annihilate(n)))
    thermal_state = gibbs(H_bath, T)
    return kronecker(bloch_state(α, β), thermal_state)
end
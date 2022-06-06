###########################
#### initial_States.jl ####
###########################

### Spins ###

## Pauli Matrices ##
σx = [[0 1];[1 0]]
σy = [[0 -im];[im 0]]
σz = [[1 0];[0 -1]]

## Spin Operators ##
sx0 = 0.5*σx
sy0 = 0.5*σy
sz0 = 0.5*σz

## Spin Coupling Operators ##
# The {θ, ϕ}[i] pairs set the direction of coupling in the ith direction #
function sc(i)
    θ = [π/2 0 π/2]
    ϕ = [0 0 π/2]
    return sx0*(sin(θ[i])*cos(ϕ[i])) + sy0*(sin(θ[i])*sin(ϕ[i])) + sz0*cos(θ[i])
end

## Bloch State ##
# This sets the initial state of the spin on the Bloch sphere, defined by the polar angle α and azimuth β #
function bloch_state()
    α = -π/2
    β = 0
    return [cos(α/2)^2 0.5*exp(-im*β)*sin(α); 0.5*exp(im*β)*sin(α) sin(α/2)^2]
end

### RC ###

## Thermal Initial State ##
function gibbs(H, T)
    n = size(H,1)
    ϵ = eigen(H).values
    P = eigen(H).vectors
    𝒵 = sum(exp(-(cfac*ϵ[i])/T) for i = 1:n)
    ρ = (1/𝒵)*Diagonal([exp(-(cfac*ϵ[i])/T) for i = 1:n])
    return P*ρ*adjoint(P)
end

### Joint Initital States ###
function ρ0(dim)
    if dim == 1
        return kronecker(bloch_state(), gibbs(HB(n1, Ω1), TDyn))
    elseif dim == 2
        return kronecker(bloch_state(), gibbs(HB(n1, Ω1), TDyn), gibbs(HB(n2, Ω2), TDyn))
    elseif dim == 3
        return kronecker(bloch_state(), gibbs(HB(n1, Ω1), TDyn), gibbs(HB(n2, Ω2), TDyn), gibbs(HB(n3, Ω3), TDyn))
    else
        print("Please return a dimension of either 1, 2 or 3.")
    end
end
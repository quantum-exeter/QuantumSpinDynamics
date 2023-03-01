#######################
#### wkCoupling.jl ####
#######################

### Weak-Coupling Integrals ###
I1(prm::LorPrm1D) = Σ(prm)
I2(prm::LorPrm1D) = 𝒬(prm) - Σ(prm)
I3(prm::LorPrm1D, β) = Δ(prm, β)

I1′(prm::LorPrm1D) = 2*Σ′(prm)
I2′(prm::LorPrm1D) = -2*Σ′(prm)
I3′(prm::LorPrm1D, β) = Δ′(prm, β)

function I4(prm::LorPrm1D)
    I(ω) = spectral_density_Lor(ω, prm)/(ω + 1)^2
    return quadgk(ω -> I(ω), 0.0, Inf)[1]
end

### Magnetisation ###
function szWK(prm::LorPrm1D, β)
    s = 1/2
    return -2*(-s*tanh(β*s) - s*(s+1)*I1′(prm) - s^2*I2′(prm) + s*tanh(β*s)*I3′(prm, β) + s^2*β*sech(β*s)^2*I3(prm, β))
end

### Zero-Temp Magnetisation ###
szWKZT(prm::LorPrm1D) = 1 - I4(prm)
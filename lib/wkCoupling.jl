#######################
#### wkCoupling.jl ####
#######################

I1(prm::LorPrm1D) = Σ(prm)
I2(prm::LorPrm1D) = 𝒬(prm) - Σ(prm)
I3(prm::LorPrm1D, β) = Δ(prm, β)

I1′(prm::LorPrm1D) = Σ′(prm)
I2′(prm::LorPrm1D) = -Σ′(prm)
I3′(prm::LorPrm1D, β) = Δ′(prm, β)

function I4′(prm::LorPrm1D)
    I(ω) = spectral_density_Lor(ω, prm)/(ω + 1)^2
    return -quadgk(ω -> I(ω), 0.0, Inf)[1]
end

𝒵(prm::LorPrm1D, β) = 2*(cosh(big(β)/2))*(1 + 0.75*big(β)*I1(prm) + 0.25*big(β)*I2(prm)) - (sinh(big(β)/2))*big(β)*I3(prm, β)
𝒵′(prm::LorPrm1D, β) = big(β)*(sinh(big(β)/2))*(1 + 0.75*big(β)*I1(prm) + 0.25*big(β)*I2(prm)) + 2*(cosh(big(β)/2))*(0.75*big(β)*I1′(prm) + 0.25*big(β)*I2′(prm)) - 0.5*big(β)^2*(cosh(big(β)/2))*I3(prm, β) - big(β)*(sinh(big(β)/2))*I3′(prm, β)

szWK(β, prm::LorPrm1D) = (1/big(β))*(1/𝒵(prm, β))*𝒵′(prm, β)
szWKZT(prm::LorPrm1D) = 1 + 2*I4′(prm)
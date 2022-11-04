#######################
#### wkCoupling.jl ####
#######################

I1(prm::LorPrm1D) = Σ(prm)
I2(prm::LorPrm1D) = 𝒬(prm) - Σ(prm)
I3(prm::LorPrm1D, β) = Δ(prm, β)

I1′(prm::LorPrm1D) = Σ′(prm)
I2′(prm::LorPrm1D) = -Σ′(prm)
I3′(prm::LorPrm1D, β) = Δ′(prm, β)

𝒵(prm::LorPrm1D, β) = 2*big(cosh(β))*(1 + 0.75*β*I1(prm) + 0.25*β*I2(prm)) - 2*big(sinh(β))*β*I3(prm, β)
𝒵′(prm::LorPrm1D, β) = 2*β*big(sinh(β))*(1 + 0.75*β*I1(prm) + 0.25*β*I2(prm)) + 2*big(cosh(β))*(3*β*I1′(prm) + β*I2′(prm)) - 2*β^2*big(cosh(β))*I3(prm, β) - 2*β*big(sinh(β))*I3′(prm, β)

szWK(β, prm::LorPrm1D) = (1/β)*(1/𝒵(prm, β))*𝒵′(prm, β)
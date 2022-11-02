#######################
#### wkCoupling.jl ####
#######################

I1(prm::LorPrm1D) = 2*Σ(prm)[1]
I2(prm::LorPrm1D) = 𝒬(prm)[1] - 2*Σ(prm)[1]
I3(prm::LorPrm1D, β) = 2*Δ(prm, β)[1]

I1′(prm::LorPrm1D) = 2*Σ′(prm)[1]
I2′(prm::LorPrm1D) = -2*Σ′(prm)[1]
I3′(prm::LorPrm1D, β) = 2*Δ′(prm, β)[1]

𝒵(prm::LorPrm1D, β) = 2*cosh(β)*(1 + 0.75*β*I1(prm) + β*I2(prm)) - 2*sinh(β)*β*I3(prm, β)
𝒵′(prm::LorPrm1D, β) = 2*β*sinh(β)*(1 + 0.75*β*I1(prm) + β*I2(prm)) + 2*cosh(β)*(0.75*β*I1′(prm) + β*I2′(prm)) - 2*β^2*cosh(β)*I3(prm, β) - 2*β*sinh(β)*I3′(prm, β)

szWK(β, prm::LorPrm1D) = (1/β)*(1/𝒵(prm, β))*𝒵′(prm, β)
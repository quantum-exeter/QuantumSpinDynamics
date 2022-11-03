#######################
#### wkCoupling.jl ####
#######################

I1(prm::LorPrm1D, ωL) = Σ(prm, ωL)
I2(prm::LorPrm1D, ωL) = 𝒬(prm) - Σ(prm, ωL)
I3(prm::LorPrm1D, β, ωL) = Δ(prm, β, ωL)

# I1′(prm::LorPrm1D) = Σ′(prm)[1]
# I2′(prm::LorPrm1D) = -Σ′(prm)[1]
# I3′(prm::LorPrm1D, β) = Δ′(prm, β)[1]

𝒵(prm::LorPrm1D, β) = 2*big(cosh(β))*(1 + 3*β*I1(prm) + β*I2(prm)) - 2*big(sinh(β))*β*I3(prm, β)
𝒵′(prm::LorPrm1D, β) = 2*β*big(sinh(β))*(1 + 3*β*I1(prm) + β*I2(prm)) + 2*big(cosh(β))*(3*β*I1′(prm) + β*I2′(prm)) - 2*β^2*big(cosh(β))*I3(prm, β) - 2*β*big(sinh(β))*I3′(prm, β)

function 𝒵(prm::LorPrm1D, β, ωL)
    2*(cosh(β*ωL/2))*(1 + 0.75*β*I1(prm, ωL) + 0.25*β*I2(prm, ωL)) - 2*(sinh(β*ωL/2))*β*I3(prm, β, ωL)
end

function szWK(β, prm::LorPrm1D)
    f(x) = log(𝒵(prm::LorPrm1D, β, x))
    df(x) = ForwardDiff.derivative(f, x)
    return (1/β)*df(1.)
end
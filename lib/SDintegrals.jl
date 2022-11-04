######################
#### integrals.jl ####
######################

function 𝒬(prm::LorPrm1D)
  return quadgk(ω -> spectral_density_Lor_inv(ω, prm), 0.0, Inf)[1]
end

function Σ(prm::LorPrm1D)
  I(ω) = spectral_density_Lor(ω, prm)*ω/(ω +1)
  Il = quadgk_cauchy(I, 0.0, 1.0, 2.0)[1]
  Ir = quadgk(ω -> I(ω)/(ω - 1), 2.0, Inf)[1]
  return Il + Ir
end
  
function Σ′(prm::LorPrm1D)
  I(ω) = spectral_density_Lor(ω, prm)*ω/(ω + 1)^2
  Il = quadgk_hadamard(I, 0.0, 1.0, 2.0)[1]
  Ir = quadgk(ω -> I(ω)/(ω - 1)^2, 2.0, Inf)[1]
  return Il + Ir
end

function Δ(prm::LorPrm1D, β)
  I(ω) = spectral_density_Lor_inv(ω, prm)*xcoth(β*ω)/(β)/(ω + 1)
  Il = quadgk_cauchy(I, 0.0, 1.0, 2.0)[1]
  Ir = quadgk(ω -> I(ω)/(ω - 1), 2.0, Inf)[1]
  return Il + Ir
end

function Δ′(prm::LorPrm1D, β)
  I(ω) = spectral_density_Lor_inv(ω, prm)*(ω^2 + 1)*xcoth(β*ω)/(β)/(ω + 1)^2
  Il = quadgk_hadamard(I, 0.0, 1.0, 2.0)[1]
  Ir = quadgk(ω -> I(ω)/(ω - 1)^2, 2.0, Inf)[1]
  return Il + Ir
end


############################
#### spectralDensity.jl ####
############################

### Spectral Density ###

spectral_density(ω, prm::LorPrm1D) = [((prm.Γ1*ω)/π)*(Λ^2/(ω^2 + Λ^2))]
spectral_density(ω, prm::LorPrm2D) = [((prm.Γ1*ω)/π)*(Λ^2/(ω^2 + Λ^2))
                                      ((prm.Γ2*ω)/π)*(Λ^2/(ω^2 + Λ^2))]
spectral_density(ω, prm::LorPrm3D) = [((prm.Γ1*ω)/π)*(Λ^2/(ω^2 + Λ^2))
                                      ((prm.Γ2*ω)/π)*(Λ^2/(ω^2 + Λ^2))
                                      ((prm.Γ3*ω)/π)*(Λ^2/(ω^2 + Λ^2))]

spectral_density_inv(ω, prm::LorPrm1D) = [((prm.Γ1)/π)*(Λ^2/(ω^2 + Λ^2))]
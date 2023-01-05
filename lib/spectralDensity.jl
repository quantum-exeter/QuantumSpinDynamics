############################
#### spectralDensity.jl ####
############################

### Spectral Density ###

## Ohmic (Algebraic) ##
spectral_density_Ohm(ω, prm::LorPrm1D) = [((prm.Γ1*ω)/π)*(Λ^2/(ω^2 + Λ^2))]
spectral_density_Ohm(ω, prm::LorPrm2D) = [((prm.Γ1*ω)/π)*(Λ^2/(ω^2 + Λ^2))
                                      ((prm.Γ2*ω)/π)*(Λ^2/(ω^2 + Λ^2))]
spectral_density_Ohm(ω, prm::LorPrm3D) = [((prm.Γ1*ω)/π)*(Λ^2/(ω^2 + Λ^2))
                                      ((prm.Γ2*ω)/π)*(Λ^2/(ω^2 + Λ^2))
                                      ((prm.Γ3*ω)/π)*(Λ^2/(ω^2 + Λ^2))]
                                      
## Lorentzian ##
spectral_density_Lor(ω, prm::LorPrm1D) = ((prm.α1*prm.Γ1*ω)/π)*(1/((prm.ω01^2 - ω^2)^2 + prm.Γ1^2*ω^2))
spectral_density_Lor_inv(ω, prm::LorPrm1D) = ((prm.α1*prm.Γ1)/π)*(1/((prm.ω01^2 - ω^2)^2 + prm.Γ1^2*ω^2))

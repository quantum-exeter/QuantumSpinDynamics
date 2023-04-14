#########################
#### entanglement.jl ####
#########################

### Entanglement Entropy ###
entropy(ρ) = -tr(ρ*log(ρ))
#########################
#### entanglement.jl ####
#########################

### Entanglement entropy ###
entropy(ρ) = -tr(ρ*log(ρ))
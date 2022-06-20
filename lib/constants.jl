######################
#### constants.jl ####
######################

γ = -1.76*10^(11) # Gyromagnetic ratio for an electron (T^-1s^-1)
Bext = 10 # External magnetic field (T)
ωL = abs(γ)*Bext # Larmor frequency (s^-1)
Λ = 10^10 # Cutoff frequency for spectral density
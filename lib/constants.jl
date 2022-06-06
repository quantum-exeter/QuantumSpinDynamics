######################
#### constants.jl ####
######################

γ = -1.76*10^(11) # Gyromagnetic ratio for an electron (T^-1s^-1)
Bext = 10 # External magnetic field (T)
ωL = abs(γ)*Bext # Larmor frequency (s^-1)
ħ = 1.05*10^(-34) # Reduced Planck constant (Js)
kB = 1.38*10^(-23) # Boltzmann cosntant (JK^-1)
Λ = 10^10 # Cutoff frequency for spectral density
cfac = (ħ*ωL)/kB
### Constants ###
γ = -1.76*10^(11) # Gyromagnetic ratio for an electron (T^-1s^-1)
Bext = 10 # External magnetic field (T)
ωL = abs(γ)*Bext # Larmor frequency (s^-1)
ħ = 1.05*10^(-34) # Reduced Planck constant (Js)
kB = 1.38*10^(-23) # Boltzmann cosntant (JK^-1)
Λ = 10^10 # Cutoff frequency for spectral density
cfac = (ħ*ωL)/kB

### Parameter Sets ###
prma = [1.4 0 10]
prmb = [7 0 50]
prmc = [0.014 0 0.1]
prmd = [1.4 0 1000]
prm = prma;

ω01, Γ1, α1 = prm*ωL
ω02, Γ2, α2 = prm*ωL
ω03, Γ3, α3 = prm*ωL

### RC-Specific Parameters ###

# Coupling Direction 1 #
Ω1 = ω01
λ1 = ωL*sqrt(α1/Ω1)
δ1 = Γ1/(2π*Ω1)
n1 = 5

# Coupling Direction 2 #
Ω2 = ω02
λ2 = ωL*sqrt(α2/Ω2)
δ2 = Γ2/(2π*Ω2)
n2 = 5

# Coupling Direction 3 #
Ω3 = ω03
λ3 = ωL*sqrt(α3/Ω3)
δ3 = Γ3/(2π*Ω3)
n3 = 5
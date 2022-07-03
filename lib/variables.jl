######################
#### variables.jl ####
######################

γ = -1 # Gyromagnetic ratio for an electron (T^-1s^-1)
Bext = 10 # External magnetic field (T)
ωL = abs(γ)*Bext # Larmor frequency (s^-1)
Λ = 10^10 # Cutoff frequency for spectral density

### Temperatures ###

### Statics ###
T = exp10.(range(-3, 3, length=100));

### Dynamics ###

TDyn = 1

### Parameter Sets ###
prma = [1.4 0.001 10]
prmb = [7 0 50]
prmc = [0.014 0 0.1]
prmd = [1.4 0 10000]
prme = [2 0.001 10]
prmf = [2 0.001 1]
prmg = [2 0.001 50]

prm1 = prme; # Change the RHS here to change parameter set for coupling direction 1
prm2 = prme; # Change the RHS here to change parameter set for coupling direction 2
prm3 = prme; # Change the RHS here to change parameter set for coupling direction 3

ω01, Γ1, α1 = prm1
ω02, Γ2, α2 = prm2
ω03, Γ3, α3 = prm3

### RC-Specific Parameters ###

Ω1, Ω2, Ω3 = [ω01 ω02 ω03]
λ1, λ2, λ3 = [sqrt(α1/Ω1) sqrt(α2/Ω2) sqrt(α3/Ω3)]
δ1, δ2, δ3 = [Γ1 Γ2 Γ3]
n1, n2, n3 = [3 3 1] # Change the number of RC levels here

δ_list(dim) = [δ1 δ2 δ3][dim] # List of dissipation strengths for Dynamics.jl module
hspace_dimension(dim) = [2*n1 2*n1*n2 2*n1*n2*n3][dim] # List of Hilbert space dimensions
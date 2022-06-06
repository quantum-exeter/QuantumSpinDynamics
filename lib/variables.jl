######################
#### variables.jl ####
######################

include("constants.jl")

### Timescales ###
ti, tf, dt = [0 10000 1000]
tspan = (ti, tf)
t = ti:dt:tf

### Temperatures ###

### Statics ###
T = exp10.(range(-3, 4, length=100));

### Dynamics ###
TDyn = 10;

### Parameter Sets ###
prma = [1.4 0 10]
prmb = [7 0 50]
prmc = [0.014 0 0.1]
prmd = [1.4 0 1000]

prm1 = prma; # Change the RHS here to change parameter set for coupling direction 1
prm2 = prma; # Change the RHS here to change parameter set for coupling direction 2
prm3 = prma; # Change the RHS here to change parameter set for coupling direction 3

ω01, Γ1, α1 = prm1*ωL
ω02, Γ2, α2 = prm2*ωL
ω03, Γ3, α3 = prm3*ωL

### RC-Specific Parameters ###

Ω1, Ω2, Ω3 = [ω01 ω02 ω03]
λ1, λ2, λ3 = ωL*[sqrt(α1/Ω1) sqrt(α2/Ω2) sqrt(α3/Ω3)]
δ1, δ2, δ3 = [Γ1/(2π*Ω1) Γ2/(2π*Ω2) Γ3/(2π*Ω3)]
n1, n2, n3 = [2 2 2] # Change the number of RC levels here

δ_list(dim) = [δ1 δ2 δ3][dim] # List of dissipation strengths for Dynamics.jl module
hspace_dimension(dim) = [2*n1 2*n1*n2 2*n1*n2*n3][dim] # List of Hilbert space dimensions
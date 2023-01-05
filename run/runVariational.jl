using CSV
using DataFrames

include("../lib/Variational.jl")
using .Variational

### Low Gamma ###
#prma = 2., 0.001, 10.
#prmb = 2., 0.001, 1.

### Weak ###
#prmc = 2., 0.001, 0.1

### High Gamma ###
#prmd = 2., 0.6, 10.
#prme = 2., 0.6, 1.

### Ultrastrong ###
#prmf = 2., 0.6, 20.
#prmg = 2., 0.6, 60.
#prmh = 2., 0.6, 80.
#prmi = 2., 0.6, 500.
#prmj = 2., 0.6, 1000.

### Parameters ###
prm = LorPrm3D(2., 0.6, 0.1, 2., 0.6, 0.1, 2., 0.6, 0.1) # Lorentzian parameters
ang =  CouplAng3D(π/2, 0.0, π/2, π/2, 0.0, 0.0) # Coupling angles
n = Lev3D(2, 2, 2) # Number of RC levels

ϕm(prm, ang, n)
ϕp(prm, ang, n)
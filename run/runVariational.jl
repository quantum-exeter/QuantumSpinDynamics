using CSV
using DataFrames

include("../lib/Statics.jl")
using .Statics

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
prm = LorPrm3D(2., 0.6, 10., 2., 0.6, 10., 2., 0.6, 10.) # Lorentzian parameters
ang =  CouplAng3D(π/2, 0.0, π/2, π/2, 0.0, 0.0) # Coupling angles
n = Lev3D(2, 2, 2) # Number of RC levels

b0 = [1, 0]
b1 = [0, 1]

using Kronecker

kronecker(b0, b1, b0)

x =ϕp(prm, ang, n)
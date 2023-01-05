using CSV
using DataFrames

include("../lib/Variational.jl")
using .Variational

### Parameters ###

## Lorentzian Spectral Density ##
# {ω_0, Γ, α}
#prma = 2., 0.001, 10.
#prmb = 2., 0.001, 1. 
#prmc = 2., 0.001, 0.1 

#prmd = 2., 0.6, 10.
#prme = 2., 0.6, 1.
#prmf = 2., 0.6, 0.1

#prmg = 2., 0.6, 20.
#prmh = 2., 0.6, 60.
#prmi = 2., 0.6, 80.
#prmj = 2., 0.6, 500.
#prmk = 2., 0.6, 1000.

prm = LorPrm3D(2., 0.6, 10., 2., 0.6, 10., 2., 0.6, 10.)

## Coupling Angles ##
ang =  CouplAng3D(π/2, 0.0, π/2, π/2, 0.0, 0.0)

## RC Levels ##
n = Lev3D(2, 2, 2) # Number of RC levels

## Temperature Range ##
T = exp10.(range(-2, 3, length=100))

ϕm(prm, ang, n)
ϕp(prm, ang, n)
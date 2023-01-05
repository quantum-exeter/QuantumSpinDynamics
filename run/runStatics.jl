using CSV
using DataFrames

include("../lib/Statics.jl")
using .Statics

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

prm = LorPrm1D(2., 0.6, 1.) 
# prm = LorPrm2D(2., 0.6, 1., 2., 0.6, 1.) 
# prm = LorPrm3D(2., 0.6, 1., 2., 0.6, 1., 2., 0.6, 1.)

## Coupling Angles ##
ang =  CouplAng1D(π/4, 0.0)
# ang =  CouplAng2D(π/2, 0.0, π/2, π/2)
# ang =  CouplAng3D(π/2, 0.0, π/2, π/2, 0.0, 0.0)

## RC Levels ##
n = Lev1D(100) # Number of RC levels
# n = Lev2D(10, 10) # Number of RC levels
# n = Lev3D(5, 5, 5) # Number of RC levels

## Temperature Range ##
T = exp10.(range(-2, 3, length=100))

sxG_list = [realIfClose(sxGibbs(i)) for i in T]
syG_list = [realIfClose(syGibbs(i)) for i in T]
szG_list = [realIfClose(szGibbs(i)) for i in T]
sxMFGS_list = [realIfClose(sxMFGS(prm, ang, n, i)) for i in T]
syMFGS_list = [realIfClose(syMFGS(prm, ang, n, i)) for i in T]
szMFGS_list = [realIfClose(szMFGS(prm, ang, n, i)) for i in T]


### Store Values ###
dfGibbs = DataFrame(hcat(T, sxG_list, syG_list, szG_list), :auto)
dfMFGS = DataFrame(hcat(T, sxMFGS_list, syMFGS_list, szMFGS_list), :auto)

### Export for Mac ###
# CSV.write("/Users/charliehogg/filename.csv",  dfMFGS, header = ["T", "sxMFGS", "syMFGS", "szMFGS"])

### Export for Windows ###
# CSV.write("C:/Users/crh222/filename.csv",  dfMFGS, header = ["T", "sxMFGS", "syMFGS", "szMFGS"])
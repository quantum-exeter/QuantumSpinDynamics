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
#prmi = 2., 0.6, 200.
#prmj = 2., 0.6, 500.

### Parameters ###
# prm = LorPrm1D(2., 0.6, 1.) # Lorentzian parameters
# prm = LorPrm2D(2.0, 0.001, 10.0, 2.0, 0.001, 10.0) # Lorentzian parameters
prm = LorPrm3D(2., 0.6, 40., 2., 0.6, 40., 2., 0.6, 40.) # Lorentzian parameters
# ang =  CouplAng1D(π/4, 0.0) # Coupling angles
# ang =  CouplAng2D(π/2, 0.0, π/2, π/2) # Coupling angles
ang =  CouplAng3D(π/2, 0.0, π/2, π/2, 0.0, 0.0) # Coupling angles
# n = Lev1D(100) # Number of RC levels
# n = Lev2D(10, 10) # Number of RC levels
n = Lev3D(5, 5, 5) # Number of RC levels

### Temperature Range ###
T = exp10.(range(-2, 3, length=100))

sxG_list = [realIfClose(sxGibbs(i)) for i in T]
syG_list = [realIfClose(syGibbs(i)) for i in T]
szG_list = [realIfClose(szGibbs(i)) for i in T]
sxMFGS_list = [realIfClose(sxMFGS(prm, ang, n, i)) for i in T]
syMFGS_list = [realIfClose(syMFGS(prm, ang, n, i)) for i in T]
szMFGS_list = [realIfClose(szMFGS(prm, ang, n, i)) for i in T]


### Store Values ###
dfGibbs = DataFrame(hcat(T, sxG_list, syG_list, szG_list), :auto)
CSV.write("Gibbs.csv",  dfGibbs, header = ["T", "sxG", "syG", "szG"])

dfMFGS = DataFrame(hcat(T, sxMFGS_list, syMFGS_list, szMFGS_list), :auto)

### Export for Mac ###
CSV.write("/Users/charliehogg/Dropbox/PhD/1. 3D Project/Data/Ultrastrong/qu_MFGS_3D_prm40_5.csv",  dfMFGS, header = ["T", "sxMFGS", "syMFGS", "szMFGS"])

### Export for Windows ###
CSV.write("C:/Users/crh222/QuantumSpinDynamics/qu_MFGS_3D_prm40_5.csv",  dfMFGS, header = ["T", "sxMFGS", "syMFGS", "szMFGS"])
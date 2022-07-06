using CSV
using DataFrames

include("../lib/Statics.jl")
using .Statics

prma = [2.0, 0.001, 10.0]

### Parameters ###
# prm = LorPrm1D(2.0, 0.001, 10.0) 
# prm = LorPrm2D(2.0, 0.001, 10.0, 2.0, 0.001, 10.0) 
prm = LorPrm3D(2.0, 0.001, 10.0, 2.0, 0.001, 10.0, 2.0, 0.001, 10.0) # Lorentzian parameters
# ang =  CouplAng1D(0.0, 0.0) # Coupling angles
# ang =  CouplAng2D(π/2, 0.0, 0.0, 0.0) # Coupling angles
ang =  CouplAng3D(π/2, 0.0, π/2, π/2, 0.0, 0.0) # Coupling angles
# n = Lev1D(100) # Number of RC levels
# n = Lev2D(10, 10) # Number of RC levels
n = Lev3D(5, 5, 5) # Number of RC levels

### Temperature Range ###
T = exp10.(range(-2, 3, length=100))

szGibbs_list = [realIfClose(szGibbs(i)) for i in T]
szMFGS_list = [realIfClose(szMFGS(prm, ang, n, i)) for i in T]

### Store Values ###
dfGibbs = DataFrame(hcat(T, szGibbs_list), :auto)
CSV.write("C://Users//crh222//Data//Quantum//Statics//Gibbs.csv",  dfGibbs, header = ["T", "szGibbs"])

dfMFGS = DataFrame(hcat(T, szMFGS_list), :auto)
CSV.write("C://Users//crh222//Data//Quantum//Statics//MFGS_3D_xyz_prma_5.csv",  dfMFGS, header = ["T", "szMFGS"])
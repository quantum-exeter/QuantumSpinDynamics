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

# MAKE SURE TO SET ALPHA AS DOUBLE THE RC VALUE #
prm = LorPrm1D(2., 0.6, 0.2)

## Temperature Range ##
T = exp10.(range(-2, 2, length=20))

szWKZT(prm)

szWK_list = [realIfClose(szWK(prm, 1/i)) for i in T]

### Store Values ###
dfWK = DataFrame(hcat(T, szWK_list), :auto)

### Export for Mac ###
CSV.write("paper_data/WK_analytical_prmc.csv",  dfWK, header = ["T", "szWK"])

### Export for Windows ###
# CSV.write("C:/Users/crh222/filename.csv",  dfMFGS, header = ["T", "szWK"])
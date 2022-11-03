using CSV
using DataFrames

include("../lib/Statics.jl")
using .Statics

#prma = 2.0, 0.001, 10.0
#prmb = 2.0, 0.001, 1.0
#prmUS = 2.0, 0.001, 1000

### Parameters ###
prm = LorPrm1D(2.0, 0.001, 1.) 

### Temperature Range ###
T = exp10.(range(-2, 3, length=100))

szWK_list = [realIfClose(szWK(1/i, prm)) for i in T]

### Store Values ###
dfWK = DataFrame(hcat(T, szWK_list), :auto)
CSV.write("C://Users//crh222//Data//3DProj//Quantum//Statics//WK.csv",  dfWK, header = ["T", "szWK"])
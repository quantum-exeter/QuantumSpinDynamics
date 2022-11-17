using CSV
using DataFrames

include("../lib/Statics.jl")
using .Statics

#prmc = 2., 0.001, 0.1
#prmf = 2., 0.6, 0.1

### Parameters ###
prm = LorPrm1D(2., 0.6, 0.1) 

### Temperature Range ###
T = exp10.(range(-3, 3, length=20))

szWK_list = [realIfClose(szWK(1/i, prm)) for i in T]

### Store Values ###
dfWK = DataFrame(hcat(T, szWK_list), :auto)
CSV.write("/Users/charliehogg/Work/QuantumSpinDynamics/paper_data/WK_prmf.csv",  dfWK, header = ["T", "szWK"])
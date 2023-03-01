using MKL
using CSV
using DataFrames
using ProgressMeter

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

# prm = LorPrm1D(2., 0.6, 0.01) 
# prm = LorPrm2D(2., 0.6, 1000., 2., 0.6, 1000.) 
prm = LorPrm3D(2., 0.6, 0.1, 2., 0.6, 0.1, 2., 0.6, 0.1)

## Coupling Angles ##
# ang =  CouplAng1D(π/4, 0.0)
# ang =  CouplAng2D(π/2, 0.0, π/2, π/2)
ang =  CouplAng3D(π/2, 0.0, π/2, π/2, 0.0, 0.0)

## RC Levels ##
# n = Lev1D(100) # Number of RC levels
# n = Lev2D(10, 10) # Number of RC levels
n = Lev3D(6, 6, 6) # Number of RC levels

## Temperature Range ##
T = exp10.(range(-2, 2, length=100))

# sxG_list = [realIfClose(sxGibbs(i)) for i in T]
# syG_list = [realIfClose(syGibbs(i)) for i in T]
# szG_list = [realIfClose(szGibbs(i)) for i in T]
# sxMFGS_list = [realIfClose(sxMFGS(prm, ang, n, i)) for i in T]
# syMFGS_list = [realIfClose(syMFGS(prm, ang, n, i)) for i in T]

szMFGS_list = zeros(length(T))
@showprogress for i in eachindex(T)
    szMFGS_list[i] = real(szMFGS(prm, ang, n, T[i]))
end

# szMFGS_list = zeros(length(T))
# @showprogress for i in eachindex(T)
#     szMFGS_list[i] = real(szMFGS(prm, ang, n, big(T[i])))
# end

### Store Values ###
# dfGibbs = DataFrame(hcat(T, szG_list), :auto)
dfMFGS = DataFrame(hcat(T, szMFGS_list), :auto)

### Export for Mac ###
CSV.write("paper_data/WK_MFGS_prmc.csv",  dfMFGS, header = ["T", "szMFGS"])

### Export for Windows ###
# CSV.write(".//paper_data//WK_MFGS_prmc.csv",  dfMFGS, header = ["T", "szMFGS"])